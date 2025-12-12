suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

## =========================
## 0) Paths & settings
## =========================
DDSDIR   <- "~/Ovary_signatures/2_1_Normalized_corrected_counts/"
OUTDIR_DE <- "~/Ovary_signatures/3_1_DE_results/"
dir.create(OUTDIR_DE, showWarnings = FALSE, recursive = TRUE)

LOG2FC_TH <- 1
PADJ_TH   <- 0.05

## =========================
## 1) Define the 4 DDS analyses
## =========================
dds_info <- list(
  
  auto = list(
    file   = "dds_tumor_CtlAuto.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_autoencoder"
  ),
  
  homolog = list(
    file   = "dds_tumor_Ctlhomolog.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_GTEx"
  ),
  
  clean_no_neg = list(
    file   = "dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_GTEx"
  ),
  
  clean_rescued = list(
    file   = "dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_GTEx"
  )
)

## =========================
## 2) Universal DE function
## =========================
run_DE <- function(tag, info) {
  message("\n========== ", tag, " ==========")
  
  dds_path <- file.path(DDSDIR, info$file)
  dds <- readRDS(dds_path)
  
  # NO volvemos a correr DESeq(); ya está corrido.
  res <- results(dds, contrast = c(info$factor, info$num, info$den))
  
  res <- res[order(res$pvalue), ]
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df %>% select(gene, everything())
  
  # DEGs
  res_DEG <- res_df %>%
    filter(!is.na(padj)) %>%
    filter(abs(log2FoldChange) >= LOG2FC_TH,
           padj < PADJ_TH)
  
  # Save
  write.table(res_df,
              file = file.path(OUTDIR_DE, paste0(tag, "_FULL.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(res_DEG,
              file = file.path(OUTDIR_DE, paste0(tag, "_DEG.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save DDS for reproducibility
  saveRDS(dds,
          file.path(OUTDIR_DE, paste0(tag, "_dds_used.rds")))
  
  invisible(res_df)
}

## =========================
## 3) Run all 4 DEAs
## =========================
for (nm in names(dds_info)) {
  run_DE(nm, dds_info[[nm]])
}
