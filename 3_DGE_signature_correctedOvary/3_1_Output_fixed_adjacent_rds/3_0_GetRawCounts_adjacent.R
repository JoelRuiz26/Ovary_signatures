###############################################
## Raw counts matrices + global PCA / t-SNE / UMAP (OVARY)
## - TCGA ovary tumors
## - Homologous ovary controls 
## - Autoencoder-selected controls (GTEx + TCGA adjacent+ TARGET)
## - CPM filtering on raw counts
## - PCA, t-SNE and UMAP using filtered raw counts
###############################################

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(edgeR)
  library(ggfortify)
  library(Rtsne)
  library(uwot)
})

options(stringsAsFactors = FALSE)

## =========================
## 0) Paths and output dir
## =========================
DIR_BASE   <- "~/Ovary_signatures/3_DGE_signature_correctedOvary/"
H5_PATH    <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"
STOP_IF_MISSING_H5 <- TRUE

if (STOP_IF_MISSING_H5) {
  stopifnot(file.exists(H5_PATH))
}

OUT_DIR <- file.path(DIR_BASE, "3_1_Output_fixed_adjacent_rds/")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

setwd(DIR_BASE)

## =========================
## 1) Load OCTAD metadata
## =========================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

## =========================
## 2) Define OVARY tumor samples (TCGA)
## =========================
case_ids_file <- file.path(DIR_BASE, "OSC_tcga_ids.RDS")
stopifnot(file.exists(case_ids_file))

muestras_ovary <- readRDS(case_ids_file)  # ~427 TCGA ovary tumors
case_ovary <- intersect(muestras_ovary, phenoDF$sample.id)
length(case_ovary)

## =========================
## 3) Homologous OVARY controls (GTEx + TCGA adjacent)
## =========================
controls_ovary <- phenoDF %>%
  dplyr::filter(
    grepl("OVARY", biopsy.site, ignore.case = TRUE),
    (data.source == "GTEX" & sample.type == "normal") |
      (data.source == "TCGA" & sample.type == "adjacent")
  ) %>%
  dplyr::pull(sample.id) %>%
  unique()

length(controls_ovary)

## =========================
## 4) Autoencoder-selected controls
## =========================
set.seed(123)  # for reproducibility
controls_autoencoder <- computeRefTissue(
  case_id      = case_ovary,
  adjacent     = TRUE,   # allow GTEx normals + TCGA/TARGET adjacents
  source       = "octad",
  control_size = 90      # 2 more because 2 ovary samples were deleted
) %>%
  as.character() %>%
  unique()

## Remove any overlap with true ovary controls
controls_autoencoder <- setdiff(controls_autoencoder, controls_ovary)
length(controls_autoencoder)  # expected ~88

## =========================
## 5) Union of all samples and load counts
## =========================
all_samples <- union(case_ovary, union(controls_ovary, controls_autoencoder))
length(all_samples)  # e.g. 602

## OCTAD counts are stored as log2(x + 1)
counts_all <- octad::loadOctadCounts(
  sample_vector = all_samples,
  type          = "counts",
  file          = H5_PATH
)

## Align metadata to the columns of the count matrix
metadata_all <- phenoDF %>%
  dplyr::filter(sample.id %in% colnames(counts_all)) %>%
  dplyr::mutate(sample.id = as.character(sample.id))

metadata_all <- metadata_all[match(colnames(counts_all), metadata_all$sample.id), , drop = FALSE]
stopifnot(identical(colnames(counts_all), metadata_all$sample.id))

## =========================
## 6) De-transform: log2(x + 1) --> x (raw counts)
## =========================
expr_raw_all <- 2^counts_all - 1
storage.mode(expr_raw_all) <- "double"
expr_raw_all <- expr_raw_all[rownames(counts_all), colnames(counts_all)]

## =========================
## 7) Group labels for downstream analyses
## =========================
metadata_all <- metadata_all %>%
  dplyr::mutate(
    group_pca = dplyr::case_when(
      sample.id %in% case_ovary           ~ "TCGA_tumor_ovary",
      sample.id %in% controls_autoencoder ~ "control_autoencoder (TCGAadjacent_gtex_target)",
      sample.id %in% controls_ovary       ~ "control_ovary_GTEx",
      TRUE                                ~ "other"
    )
  )

table(metadata_all$group_pca)

group_factor <- factor(metadata_all$group_pca)
names(group_factor) <- metadata_all$sample.id

## =========================
## 8) CPM filtering on raw counts (edgeR)
##    Keep genes with CPM > 1 in at least
##    min(group size) samples
## =========================
y <- edgeR::DGEList(counts = round(expr_raw_all), group = group_factor)
cpm_y <- edgeR::cpm(y)

min_group_size <- min(table(group_factor))
min_group_size

keep_genes_cpm <- rowSums(cpm_y > 1) >= min_group_size
sum(keep_genes_cpm)

expr_raw_all_filt <- expr_raw_all[keep_genes_cpm, , drop = FALSE]

## Remove zero-variance genes (numerical stability)
gene_var <- apply(expr_raw_all_filt, 1, var)
keep_var <- gene_var > 0

expr_raw_all_filt <- expr_raw_all_filt[keep_var, , drop = FALSE]
dim(expr_raw_all_filt)

## =========================
## 9) Subset matrices for specific contrasts (optional)
## =========================
## 9a) OVARY tumors + homologous ovary controls
samples_ovary_homolog <- union(case_ovary, controls_ovary)
samples_ovary_homolog <- intersect(samples_ovary_homolog, colnames(expr_raw_all_filt))

expr_raw_ovary_homolog <- expr_raw_all_filt[, samples_ovary_homolog, drop = FALSE]
metadata_ovary_homolog <- metadata_all %>%
  dplyr::filter(sample.id %in% samples_ovary_homolog)

## 9b) OVARY tumors + autoencoder-selected controls
samples_ovary_autoenc <- union(case_ovary, controls_autoencoder)
samples_ovary_autoenc <- intersect(samples_ovary_autoenc, colnames(expr_raw_all_filt))

expr_raw_ovary_autoenc <- expr_raw_all_filt[, samples_ovary_autoenc, drop = FALSE]
metadata_ovary_autoenc <- metadata_all %>%
  dplyr::filter(sample.id %in% samples_ovary_autoenc)

## =========================
## 10) Save matrices and metadata
## =========================
saveRDS(
  expr_raw_all_filt,
  file = file.path(OUT_DIR, "Ovary_TCGA_GTEx_AutoEncoder_rawCounts_all_filtered.rds")
)
saveRDS(
  metadata_all,
  file = file.path(OUT_DIR, "Metadata_Ovary_TCGA_GTEx_AutoEncoder_all.rds")
)

saveRDS(
  expr_raw_ovary_homolog,
  file = file.path(OUT_DIR, "Ovary_TCGA_tumor_plus_nonhomologousControls_rawCounts_filtered.rds")
)
saveRDS(
  metadata_ovary_homolog,
  file = file.path(OUT_DIR, "Metadata_Ovary_TCGA_tumor_plus_nonhomologousControls_filtered.rds")
)

saveRDS(
  expr_raw_ovary_autoenc,
  file = file.path(OUT_DIR, "Ovary_TCGA_tumor_plus_autoencoderControls_rawCounts_filtered.rds")
)
saveRDS(
  metadata_ovary_autoenc,
  file = file.path(OUT_DIR, "Metadata_Ovary_TCGA_tumor_plus_autoencoderControls_filtered.rds")
)

## =========================
## 11) Global PCA on raw filtered counts
## =========================
pca_raw <- prcomp(
  t(expr_raw_all_filt),   # samples as rows
  center = TRUE,
  scale. = TRUE
)

## Variance explained (for axis labels)
var_explained <- (pca_raw$sdev^2 / sum(pca_raw$sdev^2)) * 100

## PCA plot with ggfortify, legend on top, % variance on axes
pca_raw_plot <- autoplot(
  pca_raw,
  data   = metadata_all,
  colour = "group_pca",
  alpha  = 0.8,
  size   = 2
) +
  ggtitle("Global PCA (raw counts): OVARY tumors + controls") +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 10)
  )

pdf(file.path(OUT_DIR, "PCA_Ovary_rawCounts_autoplot.pdf"),
    width = 7, height = 6)
print(pca_raw_plot)
dev.off()

## =========================
## 12) t-SNE on PCA space (optional visualization)
## =========================
## Use first PCs from PCA on raw counts as input to t-SNE
set.seed(1234)
n_pcs_tsne <- min(50, ncol(pca_raw$x))

tsne_input <- pca_raw$x[, 1:n_pcs_tsne, drop = FALSE]

tsne_res <- Rtsne(
  tsne_input,
  dims      = 2,
  perplexity = 30,
  verbose    = FALSE,
  max_iter   = 1000
)

tsne_df <- as.data.frame(tsne_res$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$sample.id <- rownames(pca_raw$x)

tsne_df <- tsne_df %>%
  dplyr::left_join(
    metadata_all %>% dplyr::select(sample.id, group_pca),
    by = "sample.id"
  )

tsne_plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group_pca)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "t-SNE (based on PCA of raw counts): OVARY tumors + controls",
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

pdf(file.path(OUT_DIR, "tSNE_Ovary_rawCounts_PCAspace.pdf"),
    width = 7, height = 6)
print(tsne_plot)
dev.off()

## =========================
## 13) UMAP on PCA space (optional visualization)
## =========================
set.seed(5678)
n_pcs_umap <- min(50, ncol(pca_raw$x))

umap_input <- pca_raw$x[, 1:n_pcs_umap, drop = FALSE]

umap_res <- uwot::umap(
  umap_input,
  n_neighbors = 15,
  min_dist    = 0.3,
  n_components = 2,
  verbose     = FALSE
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sample.id <- rownames(pca_raw$x)

umap_df <- umap_df %>%
  dplyr::left_join(
    metadata_all %>% dplyr::select(sample.id, group_pca),
    by = "sample.id"
  )

umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = group_pca)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    title = "UMAP (based on PCA of raw counts): OVARY tumors + controls",
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

pdf(file.path(OUT_DIR, "UMAP_Ovary_rawCounts_PCAspace.pdf"),
    width = 7, height = 6)
print(umap_plot)
dev.off()


## End of script
###############################################
