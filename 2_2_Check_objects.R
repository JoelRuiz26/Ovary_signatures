suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

## =========================
## 0) Paths & dds files
## =========================
DDSDIR <- "~/Ovary_signatures/2_1_Normalized_corrected_counts/"

# Ajusta esta lista a lo que tengas en la carpeta
dds_files <- c(
  "dds_tumor_Ctlhomolog.rds",
  "dds_tumor_CtlAuto.rds",
  "dds_tumor_Ctlhomolog_CLEAN.rds",
  "dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE.rds",
  "dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives.rds"
)

## =========================
## 1) Load all dds objects
## =========================
dds_list <- list()
for (f in dds_files) {
  path_f <- file.path(DDSDIR, f)
  if (!file.exists(path_f)) {
    warning("File not found: ", path_f)
  } else {
    message("Loading: ", path_f)
    dds_list[[sub("\\.rds$", "", f)]] <- readRDS(path_f)
  }
}

## =========================
## 2) Inspect each dds
## =========================
for (nm in names(dds_list)) {
  dds <- dds_list[[nm]]
  
  cat("\n====================================\n")
  cat("Object:", nm, "\n")
  cat("Assays:", paste(assayNames(dds), collapse = ", "), "\n")
  cat("Genes:", nrow(dds), "  Samples:", ncol(dds), "\n")
  
  # Condition summary
  if ("condition" %in% colnames(colData(dds))) {
    cat("Condition table:\n")
    print(table(colData(dds)$condition))
  } else {
    cat("No 'condition' column in colData.\n")
  }
  
  # Choose which matrix to inspect: norm_clean if available, otherwise normalized counts
  if ("norm_clean" %in% assayNames(dds)) {
    mat <- assay(dds, "norm_clean")
    cat("Using assay: 'norm_clean'\n")
  } else {
    mat <- counts(dds, normalized = TRUE)
    cat("Using normalized counts from 'counts(dds, normalized=TRUE)'\n")
  }
  
  # Basic range
  rng <- range(mat, na.rm = TRUE)
  cat("Range of values:", paste(rng, collapse = "  "), "\n")
  
  # Negatives
  n_neg_values <- sum(mat < 0, na.rm = TRUE)
  n_genes_with_neg <- sum(rowSums(mat < 0, na.rm = TRUE) > 0)
  cat("Negative values:", n_neg_values,
      " (genes with ≥1 negative:", n_genes_with_neg, ")\n")
}

## Optional: keep dds_list in the environment for later use
dds_list
