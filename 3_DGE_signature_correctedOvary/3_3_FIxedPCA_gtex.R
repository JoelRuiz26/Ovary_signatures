###############################################
## PCA / t-SNE / UMAP en eigen-espacio corregido
## - Trabaja en el espacio de PCA (scores + loadings)
## - NO necesitas tocar tú la matriz cruda
## - Colorea por group_pca
###############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(Rtsne)
  library(uwot)
})

options(stringsAsFactors = FALSE)

## =========================
## 0) Paths
## =========================
DIR_BASE <- "~/Ovary_signatures/3_DGE_signature_correctedOvary/"
OUT_DIR  <- file.path(DIR_BASE, "3_1_Output_fixed_gtex_rds/")

## =========================
## 1) Cargar loadings corregidos, originales y metadata
## =========================

# Loadings originales (genes x PCs)
loadings_homolog <- readRDS(
  file.path(OUT_DIR, "Ovary_eigenspace_homolog_loadings_firstKPCs.rds")
)

# Loadings corregidos (genes x PCs)
loadings_homolog_corregidos <- readRDS(
  file.path(OUT_DIR, "Ovary_eigenspace_homolog_loadings_corrected.rds")
)

# Metadata (para colorear por group_pca)
metadata_ovary_homolog <- readRDS(
  file.path(OUT_DIR, "Metadata_Ovary_TCGA_tumor_plus_homologousControls_filtered.rds")
)

## =========================
## 2) Asegurar scores_homolog (muestras x PCs)
##    - Si el archivo existe, lo leemos
##    - Si NO existe, recalculamos el PCA homólogo y lo guardamos
## =========================

scores_file <- file.path(OUT_DIR, "Ovary_PCA_scores_homolog_rawCounts.rds")

if (file.exists(scores_file)) {
  message("Leyendo scores del PCA homólogo desde: ", scores_file)
  scores_homolog <- readRDS(scores_file)
} else {
  message("Archivo de scores no existe, recalculando PCA homólogo...")
  
  # Cargar matriz de conteos filtrados (genes x muestras: tumores + controles ovary)
  expr_raw_ovary_homolog <- readRDS(
    file.path(OUT_DIR, "Ovary_TCGA_tumor_plus_homologousControls_rawCounts_filtered.rds")
  )
  
  # PCA homólogo como en el script original (samples en filas)
  pca_homolog <- prcomp(
    t(expr_raw_ovary_homolog),
    center = TRUE,
    scale. = TRUE
  )
  
  scores_homolog <- pca_homolog$x  # muestras x PCs
  
  # Guardar para usos futuros
  saveRDS(
    scores_homolog,
    file = scores_file
  )
  
  message("Scores del PCA homólogo guardados en: ", scores_file)
}

## =========================
## 3) Alinear metadata con las muestras de scores
## =========================

metadata_ovary_homolog <- metadata_ovary_homolog %>%
  dplyr::filter(sample.id %in% rownames(scores_homolog))

metadata_ovary_homolog <- metadata_ovary_homolog[
  match(rownames(scores_homolog), metadata_ovary_homolog$sample.id),
  ,
  drop = FALSE
]

stopifnot(identical(rownames(scores_homolog),
                    metadata_ovary_homolog$sample.id))

## =========================
## 4) Transformar scores al eigen-espacio corregido
##    T = V^T %*% V_corr  (PC x PC)
##    scores_corr = scores * T
## =========================

stopifnot(
  all(dim(loadings_homolog) == dim(loadings_homolog_corregidos)),
  ncol(scores_homolog) == ncol(loadings_homolog)
)

T_mat <- t(loadings_homolog) %*% loadings_homolog_corregidos   # PCs x PCs

scores_corr <- as.matrix(scores_homolog) %*% T_mat              # muestras x PCs

## =========================
## =========================
## 5) PCA “real” del eigen-espacio corregido
##    - prcomp() sobre scores_corr (muestras x PCs)
##    - % varianza explicada en los ejes
##    - renombrar el grupo control para dejar claro que está corregido
## =========================

# PCA en el espacio corregido
pca_corr <- prcomp(
  scores_corr,   # muestras x PCs corregidos
  center = TRUE,
  scale. = TRUE
)

# Varianza explicada
var_expl_corr <- (pca_corr$sdev^2 / sum(pca_corr$sdev^2)) * 100

# Data frame con scores de PCA corregido
pca_df <- as.data.frame(pca_corr$x)
pca_df$sample.id <- rownames(pca_corr$x)

# Unir metadata
pca_df <- pca_df %>%
  dplyr::left_join(
    metadata_ovary_homolog %>% dplyr::select(sample.id, group_pca),
    by = "sample.id"
  )

# Renombrar grupo control para reflejar que es corregido SOLO en el plot
pca_df <- pca_df %>%
  dplyr::mutate(
    group_pca_plot = dplyr::recode(
      group_pca,
      "control_ovary_GTEx_only" = "control_ovary_GTEx_corrected"
    )
  )

# Plot PCA corregido con % de varianza
p_pca <- ggplot(pca_df,
                aes(x = PC1, y = PC2, color = group_pca_plot)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "PCA on corrected eigen-space: OVARY tumors + homologous controls",
    x     = paste0("PC1 (", round(var_expl_corr[1], 1), "%)"),
    y     = paste0("PC2 (", round(var_expl_corr[2], 1), "%)"),
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )

pdf(file.path(OUT_DIR, "P_2_PCA_Ovary_corrected_eigenspace_scores.pdf"),
    width = 7, height = 6)
print(p_pca)
dev.off()


## =========================
## 6) t-SNE sobre eigen-espacio corregido
## =========================
set.seed(1234)
n_pcs_tsne <- min(50, ncol(scores_corr))

tsne_input <- scores_corr[, 1:n_pcs_tsne, drop = FALSE]

tsne_res <- Rtsne(
  tsne_input,
  dims       = 2,
  perplexity = 30,
  verbose    = FALSE,
  max_iter   = 1000
)

tsne_df <- as.data.frame(tsne_res$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$sample.id <- rownames(scores_corr)

tsne_df <- tsne_df %>%
  dplyr::left_join(
    metadata_ovary_homolog %>% dplyr::select(sample.id, group_pca),
    by = "sample.id"
  )

p_tsne <- ggplot(tsne_df,
                 aes(x = tSNE1, y = tSNE2, color = group_pca)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "t-SNE (corrected eigen-space): OVARY tumors + homologous controls",
    x     = "t-SNE1",
    y     = "t-SNE2",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )

pdf(file.path(OUT_DIR, "P_2_tSNE_Ovary_corrected_eigenspace_scores.pdf"),
    width = 7, height = 6)
print(p_tsne)
dev.off()

## =========================
## 7) UMAP sobre eigen-espacio corregido
## =========================
set.seed(5678)
n_pcs_umap <- min(50, ncol(scores_corr))

umap_input <- scores_corr[, 1:n_pcs_umap, drop = FALSE]

umap_res <- uwot::umap(
  umap_input,
  n_neighbors = 15,
  min_dist    = 0.3,
  n_components = 2,
  verbose     = FALSE
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sample.id <- rownames(scores_corr)

umap_df <- umap_df %>%
  dplyr::left_join(
    metadata_ovary_homolog %>% dplyr::select(sample.id, group_pca),
    by = "sample.id"
  )

p_umap <- ggplot(umap_df,
                 aes(x = UMAP1, y = UMAP2, color = group_pca)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "UMAP (corrected eigen-space): OVARY tumors + homologous controls",
    x     = "UMAP1",
    y     = "UMAP2",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )

pdf(file.path(OUT_DIR, "P_2_UMAP_Ovary_corrected_eigenspace_scores.pdf"),
    width = 7, height = 6)
print(p_umap)
dev.off()

## End of script
###############################################
