###############################################
### DEFINE ALL PATHS HERE (único cambio)
###############################################
DIR_BASE      <- "~/Ovary_signatures/1_DGE_signature_auto/"
DIR_OUT_RDS   <- file.path(DIR_BASE, "1_1_Output_rds")
###############################################
# Load Required Libraries
###############################################
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)

setwd(DIR_BASE)

###############################################
# Download HDF5
###############################################
h5_path <- file.path("/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5")
#if (!file.exists(h5_path)) {
#  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
#  message(">> Descargando HDF5 (~3 GB) a: ", h5_path)
#  download.file(url, destfile = h5_path, mode = "wb", quiet = FALSE)
#}
#stopifnot(file.exists(h5_path))

###############################################
# Metadata (ExperimentHub: EH7274)
###############################################
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
dim(phenoDF)

###############################################
# Load samples 
###############################################
## Option 1: When you know the signature
muestras_ovary <- readRDS("0_OSC_tcga_ids.RDS")
length(muestras_ovary)
case_ovary <- intersect(muestras_ovary, phenoDF$sample.id)
if (length(case_ovary) < length(muestras_ovary)) warning("Missing some samples in OVARY")
#setdiff(muestras_ovary, case_ovary)

## Option 2: When you don't know the signature
#case_ovary <- phenoDF %>%
#  dplyr::filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
#                data.source == "TCGA",
#                sample.type != "adjacent") %>% 
#  dplyr::pull(sample.id)
#length(case_ovary)

###############################################
# Control samples (GTEx y TCGA-adjacent)
###############################################

###############################################
# Control samples (GTEx y TCGA-adjacent)
###############################################
controls_all <- computeRefTissue(
  case_id      = case_ovary,
  control_size = 88,        # TOP 88 controles
  source       = "octad",   # usa autoencoder
  adjacent     = TRUE,     # FALSE for only GTEX
  output       = FALSE      # no escribir archivos CSV
)

# Metadatos de los controles seleccionados
controls_meta <- phenoDF %>%
  dplyr::filter(sample.id %in% controls_all)
# ¿Cuáles de esos controles son de OVARY?
controls_ovary <- controls_meta %>%
  dplyr::filter(grepl("OVARY", biopsy.site, ignore.case = TRUE)) %>%
  dplyr::pull(sample.id)
table(controls_meta$biopsy.site)

cat(">> Total controls selected by autoencoder:", length(controls_all), "\n")
cat(">> Controls with biopsy.site == OVARY:", length(controls_ovary), "\n")
print(controls_ovary)

# (Opcional) Guardar controles para referencia
saveRDS(controls_all,  file.path(DIR_OUT_RDS, "1_1_controls_autoencoder_top88.rds"))

metadata <- phenoDF %>% 
  dplyr::filter(sample.id %in% c(controls_all, case_ovary))

saveRDS(metadata, file.path(DIR_OUT_RDS, "1_1_metadata_auto.rds"))
table(metadata$subtype)

###############################################
# DGE with 3 methods (solo OVARY)
###############################################
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
  
  saveRDS(
    res, 
    file.path(DIR_OUT_RDS, paste0("DE_full_", clade_label, "_", method_tag, ".rds"))
  )
  
  sig_all <- make_signature_all(res)
  sig_sig <- make_signature_sig(res, lfc_thr = lfc_thr, padj_thr = padj_thr)
  
  saveRDS(sig_all, file.path(DIR_OUT_RDS, paste0("DE_", clade_label, "_ALL_", method_tag, ".rds")))
  saveRDS(sig_sig, file.path(DIR_OUT_RDS, paste0("DE_", clade_label, "_signific_", method_tag, ".rds")))
  
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

# Ejecutar solo OVARY
for (m in c("edgeR","DESeq2","limma")) {
  run_and_save("OVARY", case_ovary, controls_all, m)
}
cat("Done. Signatures saved ALL and signific per method for OVARY\n")

###############################################
### Meta-signature
###############################################
edge   <- readRDS(file.path(DIR_OUT_RDS, "DE_OVARY_ALL_EdgeR.rds"))
deseq2 <- readRDS(file.path(DIR_OUT_RDS, "DE_OVARY_ALL_DESeq2.rds"))
limma  <- readRDS(file.path(DIR_OUT_RDS, "DE_OVARY_ALL_limma.rds"))

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
    log2FC_mean   = (log2FC_edgeR + log2FC_DESeq2 + log2FC_limma)/3,
    fisher_stat   = -2*(log(padj_edgeR) + log(padj_DESeq2) + log(padj_limma)),
    pvalue_fisher = pchisq(fisher_stat, df = 6, lower.tail = FALSE),
    padj_fisher   = p.adjust(pvalue_fisher, method = "BH"),
    regulation    = dplyr::if_else(log2FC_mean > 0, "up", "down")
  ) %>%
  dplyr::select(
    Symbol,
    log2FoldChange = log2FC_mean,
    regulation,
    padj = padj_fisher,
    log2FC_edgeR, padj_edgeR,
    log2FC_DESeq2, padj_DESeq2,
    log2FC_limma, padj_limma,
    fisher_stat, pvalue_fisher, padj_fisher
  )

saveRDS(meta_all, file.path(DIR_OUT_RDS, "Meta_DE_OVARY_ALL_auto.rds"))

meta_sig <- meta_all %>%
  dplyr::filter(
    sign(log2FC_edgeR) == sign(log2FC_DESeq2),
    sign(log2FC_edgeR) == sign(log2FC_limma),
    padj <= 0.01,
    abs(log2FoldChange) >= 1
  )

saveRDS(meta_sig, file.path(DIR_OUT_RDS, "Meta_DE_OVARY_significant_auto.rds"))

save.image("/STORAGE/csbig/jruiz/Ovary_data/Image_DiseaseSignature_OVARY_auto.RData")



#computing DE via edgeR
#loading whole octad expression data for514samples 
#OVARY | EdgeR -> TOTAL: 18431 (UP: 8946 | DOWN: 9485 ) | SIGNIFIC: 7493 (UP: 3431 | DOWN: 4062 )
#computing DE via DESeq2
#loading whole octad expression data for514samples 
#converting counts to integer mode
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#OVARY | DESeq2 -> TOTAL: 18431 (UP: 8938 | DOWN: 9493 ) | SIGNIFIC: 8552 (UP: 3915 | DOWN: 4637 )
#computing DE via limma-voom
#loading whole octad expression data for514samples 
#computing DE via limma
#OVARY | limma -> TOTAL: 18431 (UP: 9340 | DOWN: 9091 ) | SIGNIFIC: 5639 (UP: 2645 | DOWN: 2994 )
#Done. Signatures saved ALL and signific per method for OVARY




#https://doi.org/10.1214/13-AOAS683

#DOI: 10.7717/peerj.8206

