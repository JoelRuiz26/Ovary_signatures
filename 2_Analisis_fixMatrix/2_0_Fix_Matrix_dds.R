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
  condition = factor(meta_auto$group_pca)
)
rownames(colData_auto) <- meta_auto$sample.id

dds_tumor_CtlAuto <- DESeqDataSetFromMatrix(
  countData = counts_auto,
  colData   = colData_auto,
  design    = ~ condition
)

## =========================
## 3) Run DESeq and save (objetos "originales")
## =========================
dds_tumor_Ctlhomolog <- DESeq(dds_tumor_Ctlhomolog)
dds_tumor_CtlAuto    <- DESeq(dds_tumor_CtlAuto)

saveRDS(dds_tumor_Ctlhomolog, file.path(OUTDIR, "dds_tumor_Ctlhomolog.rds"))
saveRDS(dds_tumor_CtlAuto,    file.path(OUTDIR, "dds_tumor_CtlAuto.rds"))

## =========================
## 4) Normalized counts (base para la corrección)
## =========================
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
## 5) Autoencoder matrix for SVD (your original approach)
## =========================
# "X_auto": genes x samples (respetamos tu geometría)
X_auto <- auto_ctrl_mat

# Means along columns (samples)
gene_means_auto <- colMeans(X_auto)

# Centering along columns
X_auto_c <- sweep(X_auto, 2, gene_means_auto, "-")

## =========================
## 6) SVD on autoencoder controls
## =========================
sv <- svd(X_auto_c)

# Number of PCs as "noise"
k <- 2

# Loadings for first k PCs
V_k <- sv$v[, 1:k, drop = FALSE]
rownames(V_k) <- colnames(X_auto)  # sample names

## =========================
## 7) Project OVARY onto these PCs and reconstruct "noise"
## =========================
X_ovary   <- ovary_ctrl_mat
X_ovary_c <- sweep(X_ovary, 2, gene_means_auto, "-")

scores_ovary   <- X_ovary_c %*% V_k
noise_ovary_c  <- scores_ovary %*% t(V_k)

cat("Summary of noise component in OVARY controls (original geometry):\n")
print(summary(as.numeric(noise_ovary_c)))

## =========================
## 8) Clean OVARY matrix + detectar genes con negativos
## =========================
X_ovary_c_clean <- X_ovary_c - noise_ovary_c
X_ovary_clean   <- sweep(X_ovary_c_clean, 2, gene_means_auto, "+")

ovary_ctrl_mat_clean <- X_ovary_clean

# nº de muestras con valores negativos por gen
neg_counts_per_gene <- apply(ovary_ctrl_mat_clean, 1, function(x) sum(x < 0))
neg_counts_per_gene <- neg_counts_per_gene[neg_counts_per_gene > 0]

n_genes_with_neg <- length(neg_counts_per_gene)
genes_neg_names  <- names(neg_counts_per_gene)
genes_neg_ids    <- as.character(genes_neg_names)

cat("Genes con ≥1 valor negativo en controles GTEx corregidos:", n_genes_with_neg, "\n")

neg_genes_summary <- summary(
  as.numeric(ovary_ctrl_mat_clean[genes_neg_ids, , drop = FALSE])
)
neg_genes_range <- range(
  ovary_ctrl_mat_clean[genes_neg_ids, , drop = FALSE]
)

print(neg_genes_summary)
print(neg_genes_range)

## =========================
## 9) Construir matriz NORMALIZADA corregida (tumor + controles)
## =========================
norm_all_ovary_clean <- norm_all_ovary          # genes x todas las muestras
norm_all_ovary_clean[, is_ovary_ctrl] <- ovary_ctrl_mat_clean  # solo controles GTEx corregidos

## Guardar un dds "intermedio" con assay extra norm_clean (opcional)
dds_tumor_Ctlhomolog_CLEAN <- dds_tumor_Ctlhomolog
SummarizedExperiment::assay(dds_tumor_Ctlhomolog_CLEAN, "norm_clean") <- norm_all_ovary_clean

#saveRDS(dds_tumor_Ctlhomolog_CLEAN,
#        file.path(OUTDIR, "dds_tumor_Ctlhomolog_CLEAN.rds"))

## =========================
## 10) dds con CONTROLES CORREGIDOS y SIN genes negativos
##       (counts = ceil(norm_all_ovary_clean))
## =========================
keep_no_neg <- !rownames(norm_all_ovary_clean) %in% genes_neg_ids
norm_clean_no_neg <- norm_all_ovary_clean[keep_no_neg, , drop = FALSE]

# Redondear hacia arriba para usar como "counts"
counts_clean_no_neg <- ceiling(norm_clean_no_neg)

dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE <- DESeqDataSetFromMatrix(
  countData = counts_clean_no_neg,
  colData   = colData_homolog,
  design    = ~ condition
)


dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE <- estimateSizeFactors(dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE)
dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE <- estimateDispersions(dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE)

saveRDS(dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE,
        file.path(OUTDIR, "dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE.rds"))

## =========================
## 11) dds con CONTROLES CORREGIDOS + genes negativos RESCATADOS
##       (negativos con valores originales; resto corregido)
## =========================
norm_clean_rescued <- norm_all_ovary_clean
norm_clean_rescued[genes_neg_ids, ] <- norm_all_ovary[genes_neg_ids, ]

counts_clean_rescued <- ceiling(norm_clean_rescued)

dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives <- DESeqDataSetFromMatrix(
  countData = counts_clean_rescued,
  colData   = colData_homolog,
  design    = ~ condition
)


dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives <- estimateSizeFactors(dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives)
dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives <- estimateDispersions(dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives)

saveRDS(dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives,
        file.path(OUTDIR, "dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives.rds"))


