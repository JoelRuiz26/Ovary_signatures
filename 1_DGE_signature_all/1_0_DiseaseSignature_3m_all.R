# =========================================
# Load Required Libraries
# =========================================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)

setwd("~/Ovary_signatures/1_DGE_signature_all/")

# =========================================
# Download HDF5
# =========================================
h5_path <- file.path("/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5")
#if (!file.exists(h5_path)) {
#  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
#  message(">> Descargando HDF5 (~3 GB) a: ", h5_path)
#  download.file(url, destfile = h5_path, mode = "wb", quiet = FALSE)
#}
#stopifnot(file.exists(h5_path))  # ensure

# =========================================
# Metadata (ExperimentHub: EH7274)
# =========================================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
dim(phenoDF)

# =========================================
# Load samples 
# =========================================
##Option 1: When you know the signature
muestras_ovary <- readRDS("0_OSC_tcga_ids.RDS")
length(muestras_ovary) #[1] 427
# Ensure IDs existence available in OCTAD
case_ovary <- intersect(muestras_ovary, phenoDF$sample.id)
if (length(case_ovary) < length(muestras_ovary)) warning("Missing some samples in OVARY")
#setdiff(muestras_ovary, case_ovary) #[1] "TCGA-25-1870-01"

##Option 2: When you don't know the signature
#case_ovary <- phenoDF %>% 
#  dplyr::filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
#    data.source == "TCGA",
#    sample.type != "adjacent") %>% dplyr::pull(sample.id)
#length(case_ovary)
#[1] 426

# =========================================
# Control samples (GTEx y TCGA-adjacent)
# =========================================
controls_all <- phenoDF %>%
  dplyr::filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
                (data.source == "GTEX" & sample.type == "normal") |
                  (data.source == "TCGA" & sample.type == "adjacent")) %>% dplyr::pull(sample.id)
#length(controls_all) #[1] 88

metadata <- phenoDF %>% dplyr::filter(sample.id %in% c(controls_all,case_ovary))
saveRDS(metadata,
        "~/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/1_1_metadata.rds")
table(metadata$subtype)
#OVCA.Differentiated OVCA.Immunoreactive    OVCA.Mesenchymal OVCA.Proliferative 
#109                  90                 100                118 

# sample.id list by subtype
#ovca_subtypes_list <- metadata %>%
#  dplyr::filter(sample.id %in% case_ovary) %>%
#  split(.$subtype) %>%
#  lapply(function(x) x$sample.id)
#saveRDS(ovca_subtypes_list,"~/Ovary_signatures/1_DGE_signature_all/OVCA_TCGA_subtypes_ids_list.rds")

# =========================================
# DGE with 3 methods (solo OVARY)
# =========================================
lfc_thr  <- 1.0
padj_thr <- 0.01

pick_id_col <- function(df){
  if ("identifier" %in% names(df)) return("identifier")
  if ("ENSEMBL"   %in% names(df)) return("ENSEMBL")
  stop("No gene ID column found (identifier/ENSEMBL).")
}
pick_symbol_col <- function(df){
  if ("Symbol" %in% names(df)) return("Symbol")
  if ("Symbol_autho" %in% names(df)) return("Symbol_autho")
  NA_character_
}

make_signature_all <- function(res_df){
  sym_col <- pick_symbol_col(res_df)
  df <- if (!is.na(sym_col)) {
    res_df %>% dplyr::mutate(Symbol = toupper(.data[[sym_col]]))
  } else {
    res_df %>% dplyr::mutate(Symbol = .data[[pick_id_col(res_df)]])
  }
  df %>%
    
    dplyr::filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::transmute(
      Symbol,
      log2FoldChange,
      regulation = dplyr::if_else(log2FoldChange > 0, "up", "down"),
      padj
    )
}

make_signature_sig <- function(res_df, lfc_thr = 1, padj_thr = 0.01){
  make_signature_all(res_df) %>%
    dplyr::filter(abs(log2FoldChange) >= lfc_thr, padj <= padj_thr)
}

run_and_save <- function(clade_label, case_ids, ctrl_ids, method = c("edgeR","DESeq2","limma")){
  method <- match.arg(method)
  method_tag <- switch(method, edgeR="EdgeR", DESeq2="DESeq2", limma="limma")
  normalize_flag <- method %in% c("edgeR","DESeq2")
  
  res <- diffExp(
    case_id           = case_ids,
    control_id        = ctrl_ids,
    source            = "octad.whole",
    file              = h5_path,
    normalize_samples = normalize_flag,
    k                 = if (normalize_flag) 2 else 0,
    n_topGenes        = 10000,
    DE_method         = method,
    annotate          = TRUE,
    output            = FALSE
  )
  
  saveRDS(res, file.path("/home/jruiz/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/",
                         paste0("DE_full_", clade_label, "_", method_tag, ".rds")))
  
  sig_all <- make_signature_all(res)
  sig_sig <- make_signature_sig(res, lfc_thr = lfc_thr, padj_thr = padj_thr)
  
  saveRDS(sig_all, file.path("/home/jruiz/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/",
                             paste0("DE_", clade_label, "_ALL_", method_tag, ".rds")))
  saveRDS(sig_sig, file.path("/home/jruiz/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/",
                             paste0("DE_", clade_label, "_signific_", method_tag, ".rds")))
  
  cnt <- function(x) {
    tb   <- table(x)
    up   <- if ("up"   %in% names(tb)) tb[["up"]]   else 0L
    down <- if ("down" %in% names(tb)) tb[["down"]] else 0L
    c(up=as.integer(up), down=as.integer(down))
  }
  c_all <- cnt(sig_all$regulation)
  c_sig <- cnt(sig_sig$regulation)
  
  cat(clade_label, "|", method_tag,
      "-> TOTAL:",   nrow(sig_all), "(UP:", c_all["up"],   "| DOWN:", c_all["down"],   ")",
      "| SIGNIFIC:", nrow(sig_sig), "(UP:", c_sig["up"],  "| DOWN:", c_sig["down"], ")\n")
}

# Ejecutar solo para OVARY
for (m in c("edgeR","DESeq2","limma")) run_and_save("OVARY", case_ovary, controls_all, m)
cat("Done. Signatures saved ALL and signific per method for OVARY\n")



### ————— After running the 3 methods, load the three results and generate the meta-signature:

# Load the FULL (ALL) results
edge   <- readRDS("~/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/DE_OVARY_ALL_EdgeR.rds")
deseq2 <- readRDS("~/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/DE_OVARY_ALL_DESeq2.rds")
limma  <- readRDS("~/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/DE_OVARY_ALL_limma.rds")

meta_all <- edge %>%
  dplyr::select(Symbol,
                log2FC_edgeR  = log2FoldChange,
                padj_edgeR    = padj) %>%
  dplyr::inner_join(
    deseq2 %>% dplyr::select(Symbol,
                             log2FC_DESeq2 = log2FoldChange,
                             padj_DESeq2   = padj),
    by = "Symbol"
  ) %>%
  dplyr::inner_join(
    limma %>% dplyr::select(Symbol,
                            log2FC_limma  = log2FoldChange,
                            padj_limma    = padj),
    by = "Symbol"
  ) %>%
  dplyr::mutate(
    ## efecto combinado
    log2FC_mean   = (log2FC_edgeR + log2FC_DESeq2 + log2FC_limma)/3,
    ## Fisher con p ajustadas de cada método
    fisher_stat   = -2*(log(padj_edgeR) + log(padj_DESeq2) + log(padj_limma)),
    pvalue_fisher = pchisq(fisher_stat, df = 2*3, lower.tail = FALSE),
    padj_fisher   = p.adjust(pvalue_fisher, method = "BH"),
    ## NUEVO: regulation como en las tablas ALL
    regulation    = dplyr::if_else(log2FC_mean > 0, "up", "down")
  ) %>%
  ## Reordenar columnas para que se parezca a DE_*_ALL_*.rds
  dplyr::select(
    Symbol,
    log2FoldChange = log2FC_mean,  # nombre igual que en las otras ALL
    regulation,
    padj = padj_fisher,            # p combinada como padj principal
    ## columnas extra informativas
    log2FC_edgeR,   padj_edgeR,
    log2FC_DESeq2,  padj_DESeq2,
    log2FC_limma,   padj_limma,
    fisher_stat, pvalue_fisher, padj_fisher
  )

saveRDS(meta_all,
        file = "~/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/Meta_DE_OVARY_ALL.rds")

# Now the meta-significant: only genes that meet consistent direction and padj_fisher ≤ threshold
meta_sig <- meta_all %>%
  dplyr::filter(
    ## consistencia de dirección entre métodos
    sign(log2FC_edgeR)   == sign(log2FC_DESeq2),
    sign(log2FC_edgeR)   == sign(log2FC_limma),
    ## significancia según p combinada
    padj <= 0.01,
    ## efecto mínimo según log2FC combinado
    abs(log2FoldChange) >= 1
  )

saveRDS(meta_sig,
        file = "~/Ovary_signatures/1_DGE_signature_all/1_1_Output_rds/Meta_DE_OVARY_significant.rds")


save.image("/STORAGE/csbig/jruiz/Ovary_data/Image_DiseaseSignature_OVARY.RData")

#TOP_n = 10k n k=2
#computing DE via edgeR
#loading whole octad expression data for514samples 
#OVARY | EdgeR -> TOTAL: 18346 (UP: 10213 | DOWN: 8133 ) | SIGNIFIC: 9468 (UP: 5649 | DOWN: 3819 )
#computing DE via DESeq2
#loading whole octad expression data for514samples 
#converting counts to integer mode
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#OVARY | DESeq2 -> TOTAL: 18346 (UP: 10154 | DOWN: 8192 ) | SIGNIFIC: 9619 (UP: 5706 | DOWN: 3913 )
#computing DE via limma-voom
#loading whole octad expression data for514samples 
#computing DE via limma
#OVARY | limma -> TOTAL: 18346 (UP: 8688 | DOWN: 9658 ) | SIGNIFIC: 7122 (UP: 3242 | DOWN: 3880 )
#Done. Signatures saved ALL and signific per method for OVARY


#TOP_n = 15k n k=2
#computing DE via edgeR
#loading whole octad expression data for514samples 
#OVARY | EdgeR -> TOTAL: 18346 (UP: 10235 | DOWN: 8111 ) | SIGNIFIC: 9437 (UP: 5622 | DOWN: 3815 )
#computing DE via DESeq2
#loading whole octad expression data for514samples 
#converting counts to integer mode
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#OVARY | DESeq2 -> TOTAL: 18346 (UP: 10130 | DOWN: 8216 ) | SIGNIFIC: 9548 (UP: 5666 | DOWN: 3882 )
#computing DE via limma-voom
#loading whole octad expression data for514samples 
#computing DE via limma
#OVARY | limma -> TOTAL: 18346 (UP: 8688 | DOWN: 9658 ) | SIGNIFIC: 7122 (UP: 3242 | DOWN: 3880 )

#TOP_n = 12k n k=3
#computing DE via edgeR
#loading whole octad expression data for514samples 
#OVARY | EdgeR -> TOTAL: 18346 (UP: 10118 | DOWN: 8228 ) | SIGNIFIC: 8216 (UP: 4834 | DOWN: 3382 )
#computing DE via DESeq2
#loading whole octad expression data for514samples 
#converting counts to integer mode
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#OVARY | DESeq2 -> TOTAL: 18346 (UP: 10151 | DOWN: 8195 ) | SIGNIFIC: 9581 (UP: 5688 | DOWN: 3893 )
#computing DE via limma-voom
#loading whole octad expression data for514samples 
#computing DE via limma
#OVARY | limma -> TOTAL: 18346 (UP: 8688 | DOWN: 9658 ) | SIGNIFIC: 7122 (UP: 3242 | DOWN: 3880 )#
