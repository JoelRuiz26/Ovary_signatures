############################################################
## DESeq2 + VST + PCA-based noise removal (OVARY datasets)
## - 426 TCGA ovary tumors
## - 88 GTEx ovary controls (biologically matched)
## - 88 GTEx autoencoder-selected controls (heterogeneous tissues)
##
## Idea:
## 1) Build ONE DESeq2 object with 3 groups: tumor / ovary / autoencoder
## 2) VST on all samples
## 3) Use ONLY controls (ovary + autoencoder) to learn a PCA "noise subspace"
## 4) Remove the first k PCs associated with ovary vs autoencoder
## 5) Recompute PCA on: tumors + corrected GTEx ovary
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
DIR_BASE <- "~/Ovary_signatures/3_DGE_signature_correctedOvary/3_1_Output_fixed_gtex_rds/"
setwd(DIR_BASE)

expr_raw_all_int <- readRDS("3_0_1_expr_raw_autoGTEX.rds")  # genes x samples (integers)
metadata_all     <- readRDS("3_0_1_metadata_autoGTEX.rds")

stopifnot(identical(colnames(expr_raw_all_int), metadata_all$sample.id))

############################################################
## 1) Select samples and build a single DESeq2 object
##    Groups: tumor / ovary / autoencoder
############################################################

meta_sub <- metadata_all %>%
  dplyr::filter(group_pca %in% c("TCGA_tumor_ovary",
                                 "control_ovary_GTEx",
                                 "control_autoencoder_GTEx")) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      group_pca == "TCGA_tumor_ovary"        ~ "tumor",
      group_pca == "control_ovary_GTEx"      ~ "ovary",
      group_pca == "control_autoencoder_GTEx"~ "autoencoder",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group))

counts_sub <- expr_raw_all_int[, meta_sub$sample.id, drop = FALSE]

# Basic sanity
stopifnot(identical(colnames(counts_sub), meta_sub$sample.id))

## Global unified gene filtering (same genes for everything)
n_samples_min_global <- min(table(meta_sub$group))  # min over tumor/ovary/autoencoder
keep_global <- rowSums(counts_sub >= 10) >= n_samples_min_global
counts_sub  <- counts_sub[keep_global, , drop = FALSE]

message("Genes kept after unified filtering: ", nrow(counts_sub))

# Build DESeq2 object with 3 groups
colData_all <- data.frame(
  sample_id = meta_sub$sample.id,
  group     = factor(meta_sub$group,
                     levels = c("ovary", "autoencoder", "tumor"))
)
rownames(colData_all) <- colData_all$sample_id
colData_all$sample_id <- NULL

dds_all <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = colData_all,
  design    = ~ group
)

## Estimate size factors + VST
dds_all <- estimateSizeFactors(dds_all)
vsd_all <- vst(dds_all, blind = FALSE)

mat_vst <- assay(vsd_all)  # genes x samples

############################################################
## 2) Define indexes for each group
############################################################

group_vec <- colData(vsd_all)$group
idx_tumor <- which(group_vec == "tumor")
idx_ovary <- which(group_vec == "ovary")
idx_auto  <- which(group_vec == "autoencoder")

message("Tumor samples: ", length(idx_tumor))
message("GTEx ovary    : ", length(idx_ovary))
message("GTEx autoenc. : ", length(idx_auto))

############################################################
## 3) PCA on controls only (ovary + autoencoder)
##    This PCA defines the "noise subspace"
############################################################

# Controls-only VST matrix in samples x genes
X_ctrl <- t(mat_vst[, c(idx_ovary, idx_auto), drop = FALSE])
group_ctrl <- c(rep("ovary", length(idx_ovary)),
                rep("autoencoder", length(idx_auto)))

# PCA (samples x genes)
pca_ctrl <- prcomp(X_ctrl, center = TRUE, scale. = FALSE)

# Quick check of variance explained (you can print or plot if you want)
var_explained <- pca_ctrl$sdev^2 / sum(pca_ctrl$sdev^2)
message("PC1 variance (controls): ", round(100 * var_explained[1], 1), "%")
message("PC2 variance (controls): ", round(100 * var_explained[2], 1), "%")

# OPTIONAL: inspect separation of ovary vs autoencoder by PC1/PC2
# (commented out to keep script clean)
# df_ctrl_pca <- data.frame(PC1 = pca_ctrl$x[,1],
#                           PC2 = pca_ctrl$x[,2],
#                           group = group_ctrl)
# ggplot(df_ctrl_pca, aes(PC1, PC2, color = group)) +
#   geom_point() + theme_bw()

############################################################
## 4) Remove k PCs interpreted as "noise"
##    (e.g. PCs dominated by ovary vs autoencoder differences)
############################################################

k_noise <- 2  # <- YOU choose: 1, 2, 3, ... (inspect PCA to justify)

# Components to use as noise
loadings_noise <- pca_ctrl$rotation[, 1:k_noise, drop = FALSE]  # genes x k
center_vec     <- pca_ctrl$center                               # length = genes

## Helper: function to remove those PCs from any sample x gene matrix
remove_noise_pcs <- function(X, center_vec, loadings_noise) {
  # X: samples x genes, same genes/order as used in pca_ctrl
  X_centered <- sweep(X, 2, center_vec, FUN = "-")
  scores     <- X_centered %*% loadings_noise              # samples x k
  recon      <- scores %*% t(loadings_noise)               # samples x genes
  X_corr_cent <- X_centered - recon
  X_corrected <- sweep(X_corr_cent, 2, center_vec, FUN = "+")
  return(X_corrected)
}

############################################################
## 5) Build matrix (tumor + ovary) and remove noise PCs
############################################################

# Matrix to correct: tumor + ovary only (samples x genes)
X_tumor_ovary <- t(mat_vst[, c(idx_tumor, idx_ovary), drop = FALSE])

# Make sure gene order matches pca_ctrl rotation
stopifnot(identical(colnames(X_tumor_ovary), rownames(loadings_noise)))

# Remove noise PCs
X_tumor_ovary_corr <- remove_noise_pcs(
  X          = X_tumor_ovary,
  center_vec = center_vec,
  loadings_noise = loadings_noise
)

# Back to genes x samples
mat_corr <- t(X_tumor_ovary_corr)
colnames(mat_corr) <- colnames(mat_vst)[c(idx_tumor, idx_ovary)]
rownames(mat_corr) <- rownames(mat_vst)

############################################################
## 6) PCA on corrected matrix (tumor + corrected GTEx ovary)
############################################################

# prcomp expects samples x genes
pca_final <- prcomp(t(mat_corr), center = TRUE, scale. = FALSE)

pc_var <- pca_final$sdev^2 / sum(pca_final$sdev^2)

df_pca_final <- data.frame(
  PC1  = pca_final$x[, 1],
  PC2  = pca_final$x[, 2],
  group = factor(c(rep("tumor", length(idx_tumor)),
                   rep("ovary_corrected", length(idx_ovary)))
  )
)

p_pca_final <- ggplot(df_pca_final,
                      aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  xlab(paste0("PC1: ", round(100 * pc_var[1], 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * pc_var[2], 1), "% variance")) +
  ggtitle(paste0("PCA: TCGA tumors vs GTEx ovary (", 
                 k_noise, " PCs removed from control space)"))

ggsave("3_2_1_PCA_vst_Tumor_vs_GTEx_ovary_corrected_by_autoencoder.pdf",
       p_pca_final, width = 8, height = 6, dpi = 300)

############################################################
## 7) Save objects for downstream analyses
############################################################

saveRDS(dds_all,   file = "3_2_0_dds_all_Tumor_GTEx_ovary_autoencoder.rds")
saveRDS(vsd_all,   file = "3_2_0_vsd_all_Tumor_GTEx_ovary_autoencoder.rds")
saveRDS(mat_corr,  file = "3_2_1_vst_Tumor_plus_GTEx_ovary_corrected_by_autoencoder.rds")

message(">> Finished: unified DESeq2, control-PC noise removal, and corrected PCA.")
