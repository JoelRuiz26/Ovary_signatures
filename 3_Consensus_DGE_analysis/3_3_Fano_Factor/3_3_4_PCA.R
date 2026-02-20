suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

OUT_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor"

# -------------------------
# Load normalized count matrices + metadata
# -------------------------
norm_all_ovary <- readRDS(file.path(OUT_DIR, "3_3_2_norm_counts_tumor_GTEX.rds"))  # tumor + GTEx/adj controls
norm_all_auto  <- readRDS(file.path(OUT_DIR, "3_3_2_norm_count_tumor_AE.rds"))     # tumor + AE controls
meta_all       <- readRDS(file.path(OUT_DIR, "3_3_1_metadata_all.rds"))

# -------------------------
# PCA helper (log2(x+1) on normalized counts)
# -------------------------
do_pca_df <- function(mat_norm, meta_all, contrast_label, keep_groups) {
  # keep only samples present in this matrix, in the same order
  meta <- meta_all %>%
    filter(sample.id %in% colnames(mat_norm)) %>%
    mutate(sample.id = as.character(sample.id)) %>%
    dplyr::slice(match(colnames(mat_norm), sample.id))
  stopifnot(identical(meta$sample.id, colnames(mat_norm)))
  
  # keep ONLY the 2 groups for this contrast (this makes PCA separate per contrast)
  keep_idx <- meta$group_pca %in% keep_groups
  meta <- meta[keep_idx, , drop = FALSE]
  mat_norm <- mat_norm[, meta$sample.id, drop = FALSE]
  
  # log-transform then PCA on samples
  X <- log2(mat_norm + 1)
  pca <- prcomp(t(X), center = TRUE, scale. = FALSE)
  
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  
  tibble(
    sample.id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group = meta$group_pca,
    contrast = contrast_label,
    PC1_var = var_exp[1],
    PC2_var = var_exp[2]
  )
}

df_pca <- bind_rows(
  do_pca_df(
    norm_all_ovary, meta_all,
    "Tumor vs GTEx controls",
    keep_groups = c("TCGA_tumor", "control_GTEx")
  ),
  do_pca_df(
    norm_all_auto, meta_all,
    "Tumor vs AE controls",
    keep_groups = c("TCGA_tumor", "control_autoencoder")
  )
)

# -------------------------
# Grid PCA plot
# -------------------------
p_pca <- ggplot(df_pca, aes(PC1, PC2, color = group)) +
  geom_point(size = 2, alpha = 0.85) +
  facet_wrap(~ contrast, ncol = 2) +
  labs(
    title = "PCA of normalized counts",
    x = "PC1",
    y = "PC2",
    color = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# add per-panel % variance (computed per contrast)
df_var <- df_pca %>%
  group_by(contrast) %>%
  summarise(
    PC1 = dplyr::first(PC1_var),
    PC2 = dplyr::first(PC2_var),
    .groups = "drop"
  ) %>%
  mutate(lbl = sprintf("PC1=%.1f%% | PC2=%.1f%%", 100*PC1, 100*PC2))

p_pca +
  geom_label(
    data = df_var,
    aes(x = -Inf, y = Inf, label = lbl),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 4,
    label.size = 0
  )

