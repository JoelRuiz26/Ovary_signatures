# =========================================
# Load Required Libraries
# =========================================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)
library(harmonicmeanp)  


setwd("~/Ovary_signatures/1_DGE_signature_ovary/")

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
        "1_1_Output_rds/1_1_metadata.rds")
table(metadata$subtype)
#OVCA.Differentiated OVCA.Immunoreactive    OVCA.Mesenchymal OVCA.Proliferative 
#109                  90                 100                118 

# sample.id list by subtype
#ovca_subtypes_list <- metadata %>%
#  dplyr::filter(sample.id %in% case_ovary) %>%
#  split(.$subtype) %>%
#  lapply(function(x) x$sample.id)
#saveRDS(ovca_subtypes_list,"OVCA_TCGA_subtypes_ids_list.rds")

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
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      !is.na(pvalue)          # (OPCIONAL pero recomendable)
    ) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::transmute(
      Symbol,
      log2FoldChange,
      regulation = dplyr::if_else(log2FoldChange > 0, "up", "down"),
      pvalue,               
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
  
  saveRDS(res, file.path("1_1_Output_rds/",
                         paste0("DE_full_", clade_label, "_", method_tag, ".rds")))
  
  sig_all <- make_signature_all(res)
  sig_sig <- make_signature_sig(res, lfc_thr = lfc_thr, padj_thr = padj_thr)
  
  saveRDS(sig_all, file.path("1_1_Output_rds/",
                             paste0("DE_", clade_label, "_ALL_", method_tag, ".rds")))
  saveRDS(sig_sig, file.path("1_1_Output_rds/",
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
### 1. Load results
edge   <- readRDS("1_1_Output_rds/DE_OVARY_ALL_EdgeR.rds")
deseq2 <- readRDS("1_1_Output_rds/DE_OVARY_ALL_DESeq2.rds")
limma  <- readRDS("1_1_Output_rds/DE_OVARY_ALL_limma.rds")

### 2. Merge tables (logFC + p-values) 
# https://www.pnas.org/doi/full/10.1073/pnas.1814092116
# Para combinar la evidencia de los tres métodos de expresión diferencial
# (edgeR, DESeq2 y limma), utilizamos el método inverse normal / Stouffer.
# Primero convertimos los p-values de cada método en Z-scores con signo,
# usando el signo del logFC. Después combinamos los Z por gen:
#   Z_meta = sum(Z_i) / sqrt(k)
# y obtenemos un meta p-value bilateral:
#   p_meta = 2 * pnorm(-abs(Z_meta)).
# Finalmente, corregimos por múltiples pruebas con Benjamini–Hochberg.

meta_tmp <- edge %>%
  dplyr::select(
    Symbol,
    log2FC_edgeR = log2FoldChange,
    p_edgeR      = pvalue,
    padj_edgeR   = padj
  ) %>%
  dplyr::inner_join(
    deseq2 %>% dplyr::select(
      Symbol,
      log2FC_DESeq2 = log2FoldChange,
      p_DESeq2      = pvalue,
      padj_DESeq2   = padj
    ),
    by = "Symbol"
  ) %>%
  dplyr::inner_join(
    limma %>% dplyr::select(
      Symbol,
      log2FC_limma = log2FoldChange,
      p_limma      = pvalue,
      padj_limma   = padj
    ),
    by = "Symbol"
  )

### 3. Construir matrices de p-values y logFC (genes x métodos)
pmat <- as.matrix(meta_tmp[, c("p_edgeR", "p_DESeq2", "p_limma")])
lfc_mat <- as.matrix(meta_tmp[, c("log2FC_edgeR", "log2FC_DESeq2", "log2FC_limma")])

# Evitar problemas numéricos con p extremadamente pequeños o cercanos a 1
eps <- .Machine$double.xmin
pmat[pmat < eps]       <- eps
pmat[pmat > 1 - eps]   <- 1 - eps

### 4. Calcular Z-scores con signo para cada método
# p de entrada se asumen bilaterales; construimos Z firmado:
# Z_i = sign(logFC_i) * qnorm(p_i/2, lower.tail = FALSE)
z_mat <- sign(lfc_mat) * qnorm(pmat / 2, lower.tail = FALSE)

### 5. Z combinado por gen (Stouffer) y meta p-value bilateral
k_methods <- ncol(z_mat)

z_meta <- rowSums(z_mat, na.rm = TRUE) / sqrt(k_methods)
p_meta <- 2 * pnorm(-abs(z_meta))  # bilateral
p_meta[p_meta < eps] <- eps


### 6. Construir meta_all con Stouffer
meta_all <- meta_tmp %>%
  dplyr::mutate(
    z_meta       = z_meta,
    pvalue_meta  = p_meta,
    padj_meta    = p.adjust(pvalue_meta, method = "BH"),
    log2FC_mean  = (log2FC_edgeR + log2FC_DESeq2 + log2FC_limma) / 3,
    regulation   = dplyr::if_else(log2FC_mean > 0, "up", "down"),
    # cuántos métodos individuales son significativos (FDR < 0.01)
    n_sig_methods = (padj_edgeR  < 0.01) +
      (padj_DESeq2 < 0.01) +
      (padj_limma  < 0.01)
  ) %>%
  dplyr::select(
    Symbol,
    log2FoldChange = log2FC_mean,
    regulation,
    z_meta,
    pvalue_meta,
    padj = padj_meta,
    n_sig_methods,
    log2FC_edgeR,  p_edgeR,  padj_edgeR,
    log2FC_DESeq2, p_DESeq2, padj_DESeq2,
    log2FC_limma,  p_limma,  padj_limma
  )
meta_all$padj[meta_all$padj < eps] <- eps


saveRDS(
  meta_all,
  "1_1_Output_rds/1_2_Meta_DE_OVARY_ALL_Stouffer.rds"
)

meta_all_compare <- meta_all %>% dplyr::select(Symbol,pvalue_meta,padj,n_sig_methods, p_edgeR,p_DESeq2,p_limma)

### 7. Meta-significant & consistent genes (consenso + efecto)
meta_sig <- meta_all %>%
  dplyr::filter(
    # consistencia de dirección entre métodos
    sign(log2FC_edgeR)   == sign(log2FC_DESeq2),
    sign(log2FC_edgeR)   == sign(log2FC_limma),
    # significancia según p combinada (Stouffer)
    padj <= 0.01,
    # al menos 2 métodos individuales significativos (puedes cambiar a ==3 si quieres ultra estricto)
    n_sig_methods >= 2,
    # tamaño de efecto mínimo
    abs(log2FoldChange) >= 1
  )

saveRDS(
  meta_sig,
  "1_1_Output_rds/1_2_Meta_DE_OVARY_significant_Stouffer.rds"
)


save.image("/STORAGE/csbig/jruiz/Ovary_data/1_Output_RawOvary/Image_DiseaseSignature_OVARY.RData")
#load("/STORAGE/csbig/jruiz/Ovary_data/Image_DiseaseSignature_OVARY.RData")

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
