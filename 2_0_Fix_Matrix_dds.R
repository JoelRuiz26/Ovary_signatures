suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(SummarizedExperiment)
})

## =========================
## 0) Paths & data
## =========================
INPUTDIR <- "~/Ovary_signatures/1_1_Counts_Raw/"
OUTDIR   <- "~/Ovary_signatures/2_1_Normalized_corrected_counts/"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

expr_raw_Tumor_CtlAuto <- readRDS(file.path(INPUTDIR, "expr_raw_Tumor_CtlAuto.rds"))
counts_homolog_raw     <- readRDS(file.path(INPUTDIR, "expr_raw_Tumor_Ctlhomolog.rds"))
metadata_all           <- readRDS(file.path(INPUTDIR, "metadata_all.rds"))

## =========================
## 1) Tumor + homolog controls (dds_tumor_Ctlhomolog)
## =========================
meta_homolog <- metadata_all %>%
  dplyr::filter(sample.id %in% colnames(counts_homolog_raw))

# Reorder metadata to match count columns
meta_homolog <- meta_homolog[match(colnames(counts_homolog_raw),
                                   meta_homolog$sample.id),
                             , drop = FALSE]

# Gene filter: >= 10 counts in at least min(#samples per group)
n_min_homolog <- min(table(meta_homolog$group_pca))
keep_homolog  <- rowSums(counts_homolog_raw >= 10) >= n_min_homolog
counts_homolog <- counts_homolog_raw[keep_homolog, , drop = FALSE]

# colData
colData_homolog <- data.frame(
  condition = factor(meta_homolog$group_pca)
)
rownames(colData_homolog) <- meta_homolog$sample.id

dds_tumor_Ctlhomolog <- DESeqDataSetFromMatrix(
  countData = counts_homolog,
  colData   = colData_homolog,
  design    = ~ condition
)

## =========================
## 2) Tumor + autoencoder controls (dds_tumor_CtlAuto)
## =========================
meta_auto <- metadata_all %>%
  dplyr::filter(sample.id %in% colnames(expr_raw_Tumor_CtlAuto))

meta_auto <- meta_auto[match(colnames(expr_raw_Tumor_CtlAuto),
                             meta_auto$sample.id),
                       , drop = FALSE]

common_genes <- rownames(counts_homolog)

counts_auto <- expr_raw_Tumor_CtlAuto[common_genes,
                                      meta_auto$sample.id,
                                      drop = FALSE]

colData_auto <- data.frame(
  condition = factor(meta_auto$group_pca))

rownames(colData_auto) <- meta_auto$sample.id

dds_tumor_CtlAuto <- DESeqDataSetFromMatrix(
  countData = counts_auto,
  colData   = colData_auto,
  design    = ~ condition)

## =========================
## 3) Run DESeq and save
## =========================
dds_tumor_Ctlhomolog <- DESeq(dds_tumor_Ctlhomolog)
dds_tumor_CtlAuto    <- DESeq(dds_tumor_CtlAuto)

saveRDS(dds_tumor_Ctlhomolog, file.path(OUTDIR, "dds_tumor_Ctlhomolog.rds"))
saveRDS(dds_tumor_CtlAuto,    file.path(OUTDIR, "dds_tumor_CtlAuto.rds"))

## =========================
## Extract normalized counts
## =========================
# Matrices genes x samples
norm_all_ovary <- counts(dds_tumor_Ctlhomolog, normalized = TRUE)
norm_all_auto  <- counts(dds_tumor_CtlAuto,    normalized = TRUE)

# Condition flags
is_ovary_ctrl       <- colData(dds_tumor_Ctlhomolog)$condition == "control_GTEx"
is_autoencoder_ctrl <- colData(dds_tumor_CtlAuto)$condition    == "control_autoencoder"

# Subsets: controls only (genes x samples)
ovary_ctrl_mat <- as.matrix(norm_all_ovary[, is_ovary_ctrl,       drop = FALSE])  # GTEx ovary / homolog
auto_ctrl_mat  <- as.matrix(norm_all_auto [, is_autoencoder_ctrl, drop = FALSE])  # autoencoder

# Check that genes match
stopifnot(identical(rownames(ovary_ctrl_mat), rownames(auto_ctrl_mat)))

## =========================
## 4) Autoencoder matrix for SVD (your original approach)
## =========================
# NOTE: we keep your exact geometry: matrices treated as samples x genes
# even though the object is genes x samples. We do NOT transpose.

# "X_auto": you labeled this as samples x genes, we keep it as-is
X_auto <- auto_ctrl_mat   # actually genes x samples, but we respect your pipeline

# Means along columns (in your original code)
gene_means_auto <- colMeans(X_auto)  # length = number of autoencoder samples

# Centering along columns (your original convention)
X_auto_c <- sweep(X_auto, 2, gene_means_auto, "-")

## =========================
## 5) SVD on autoencoder controls (as in your code)
## =========================
sv <- svd(X_auto_c)

# Number of components to treat as "noise"
k <- 2

# Loadings for the first k PCs (you originally used sv$v)
V_k <- sv$v[, 1:k, drop = FALSE]
rownames(V_k) <- colnames(X_auto)  # sample names

## =========================
## 6) Project OVARY onto these PCs and reconstruct "noise"
## =========================
# OVARY matrix (you also treat this as samples x genes)
X_ovary <- ovary_ctrl_mat  # same orientation as X_auto

# Center OVARY using the same means (your original code)
X_ovary_c <- sweep(X_ovary, 2, gene_means_auto, "-")

# Scores of OVARY in the subspace of the first k PCs
scores_ovary <- X_ovary_c %*% V_k   # dim: genes x k (assuming same #controls)

# Part of X_ovary_c explained only by these k PCs
noise_ovary_c <- scores_ovary %*% t(V_k)   # dim: genes x samples

cat("Summary of noise component in OVARY controls (original geometry):\n")
print(summary(as.numeric(noise_ovary_c)))

## =========================
## 7) Get "clean" OVARY matrix (your logic)
## =========================
# Remove the part explained by the k PCs
X_ovary_c_clean <- X_ovary_c - noise_ovary_c

# Back to original scale adding the same column means
X_ovary_clean <- sweep(X_ovary_c_clean, 2, gene_means_auto, "+")

# This is your cleaned OVARY control matrix (genes x samples)
ovary_ctrl_mat_clean <- X_ovary_clean

## =========================
## 8) Count genes with negative values (fixed logic only)
## =========================
# For each gene (row), count how many samples have negative values
neg_counts_per_gene <- apply(ovary_ctrl_mat_clean, 1, function(x) sum(x < 0))

# Keep only genes with at least one negative value
neg_counts_per_gene <- neg_counts_per_gene[neg_counts_per_gene > 0]

n_genes_with_neg <- length(neg_counts_per_gene)
genes_neg_names  <- names(neg_counts_per_gene)


# Number of genes and their IDs (as character)
n_genes_with_neg <- length(neg_counts_per_gene)
genes_neg_names  <- names(neg_counts_per_gene)
genes_neg_ids    <- as.character(genes_neg_names)   # explicit character vector

# Inspect how many genes have ≥1 negative
n_genes_with_neg
head(genes_neg_ids)
sort(neg_counts_per_gene)

# Inspect dynamic range / summary of counts for these genes
neg_genes_summary <- summary(
  as.numeric(ovary_ctrl_mat_clean[genes_neg_ids, , drop = FALSE])
)
neg_genes_range <- range(
  ovary_ctrl_mat_clean[genes_neg_ids, , drop = FALSE]
)

neg_genes_summary  # min, quartiles, median, max over all values of negative genes
neg_genes_range  

## =========================
## 9) Put cleaned OVARY controls back into full matrix
##    and save a CLEAN DESeq2 object
## =========================
# Start from the full normalized matrix (tumor + controls)
norm_all_ovary_clean <- norm_all_ovary   # genes x all_samples

# Replace only the OVARY control columns with the cleaned values
norm_all_ovary_clean[, is_ovary_ctrl] <- ovary_ctrl_mat_clean

# Copy DESeq2 object and attach the cleaned norm matrix as an extra assay
dds_tumor_Ctlhomolog_CLEAN <- dds_tumor_Ctlhomolog
SummarizedExperiment::assay(dds_tumor_Ctlhomolog_CLEAN, "norm_clean") <- norm_all_ovary_clean

saveRDS(dds_tumor_Ctlhomolog_CLEAN, file.path(OUTDIR, "dds_tumor_Ctlhomolog_CLEAN.rds"))

# dds sin genes negativos
saveRDS(dds_tumor_Ctlhomolog_CLEAN[!rownames(dds_tumor_Ctlhomolog_CLEAN) %in% genes_neg_ids, ],
        file.path(OUTDIR, "dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE.rds"))

# sustituir en norm_clean los genes negativos por los valores originales
norm_clean_fix <- SummarizedExperiment::assay(dds_tumor_Ctlhomolog_CLEAN, "norm_clean")
norm_clean_fix[genes_neg_ids, ] <- norm_all_ovary[genes_neg_ids, ]
SummarizedExperiment::assay(dds_tumor_Ctlhomolog_CLEAN, "norm_clean") <- norm_clean_fix

# dds con todos los genes, negativos rescatados con valores originales
saveRDS(dds_tumor_Ctlhomolog_CLEAN,
        file.path(OUTDIR, "dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives.rds"))


