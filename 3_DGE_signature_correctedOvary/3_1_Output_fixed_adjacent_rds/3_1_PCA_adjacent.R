############################################################
## DESeq2 objects + VST-PCA for OVARY datasets
## - Matrix 1: TCGA tumors + GTEx ovary controls 
## - Matrix 2: TCGA tumors + GTEx autoencoder controls
## - Only normalization + VST + PCA 
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

## =========================
## 0) Load inputs
## =========================
DIR_BASE <- "~/Ovary_signatures/3_DGE_signature_correctedOvary/3_1_Output_fixed_adjacent_rds/"
setwd(DIR_BASE)

expr_raw_all_int <- readRDS("3_0_1_expr_raw_autoADJAC.rds")
metadata_all     <- readRDS("3_0_1_metadata_autoADJAC.rds")
table(metadata_all$group_pca)
stopifnot(identical(colnames(expr_raw_all_int), metadata_all$sample.id))

############################################################
## Function to build a DESeq2 VST-PCA object for a group pair
############################################################
process_group <- function(meta, counts, condition_levels,
                          dds_file, vsd_file, pca_file) {
  
  # Build colData
  sampleTable <- data.frame(
    condition = factor(meta$condition, levels = condition_levels)
  )
  rownames(sampleTable) <- meta$sample.id
  
  # DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = sampleTable,
    design    = ~ condition
  )
  
  # Filtering: keep genes with >=10 counts in at least
  # the minimum number of samples across conditions
  n_samples_min <- min(table(sampleTable$condition))
  keep <- rowSums(counts(dds) >= 10) >= n_samples_min
  dds  <- dds[keep, ]
  
  # Normalization + variance-stabilizing transform
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = FALSE)
  
  # PCA
  pca_plot <- plotPCA(vsd, intgroup = "condition")
  ggsave(pca_file, pca_plot, width = 8, height = 6, dpi = 300)
  
  # Save objects
  saveRDS(dds, dds_file)
  saveRDS(vsd, vsd_file)
}

############################################################
## 1) TCGA tumors + GTEx OVARY controls
############################################################

meta_ovary <- metadata_all %>%
  filter(group_pca %in% c("TCGA_tumor_ovary", "control_ovary_GTEx")) %>%
  mutate(
    condition = ifelse(group_pca == "TCGA_tumor_ovary",
                       "tumor", "control_ovary")
  )

counts_ovary <- expr_raw_all_int[, meta_ovary$sample.id, drop = FALSE]

process_group(
  meta          = meta_ovary,
  counts        = counts_ovary,
  condition_levels = c("control_ovary", "tumor"),
  dds_file      = "3_1_1_dds_Tumor_vs_GTEx_ovary.rds",
  vsd_file      = "3_1_2_vsd_Tumor_vs_GTEx_ovary.rds",
  pca_file      = "3_1_2_PCA_vst_Tumor_vs_GTEx_ovary.pdf"
)

############################################################
## 2) TCGA tumors + GTEx AUTOENCODER controls
############################################################

meta_auto <- metadata_all %>%
  filter(group_pca %in% c("TCGA_tumor_ovary", "control_autoencoder_adjacent")) %>%
  mutate(
    condition = ifelse(group_pca == "TCGA_tumor_ovary",
                       "tumor", "control_autoencoder")
  )

counts_auto <- expr_raw_all_int[, meta_auto$sample.id, drop = FALSE]

process_group(
  meta          = meta_auto,
  counts        = counts_auto,
  condition_levels = c("control_autoencoder", "tumor"),
  dds_file      = "3_1_3_dds_Tumor_vs_adjacent_autoencoder.rds",
  vsd_file      = "3_1_4_vsd_Tumor_vs_adjacent_autoencoder.rds",
  pca_file      = "3_1_4_PCA_vst_Tumor_vs_adjacent_autoencoder.pdf"
)

message(">> Finished: DESeq2 objects and VST-PCA for both datasets.")
