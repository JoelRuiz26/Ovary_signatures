# =============================================================================
# ovary_mra.R
# Master Regulator Analysis (MRA) for ovarian cancer
#
# Overview:
#   1. Cross-source concordance: correlate RRA (GEO meta-analysis) vs DESeq2
#      Wald statistics (TCGA vs AE; TCGA vs GTEx)
#   2. Core consensus signature construction: retain genes with concordant
#      direction across all three sources; apply elbow cutoff on |z_mean|
#   3. Regulon construction from ARACNe-AP network + TCGA expression matrix
#   4. MRA via msVIPER with signature-label permutation null model
#      (core signature + independent AE-z and GTEx-z signatures)
#   5. Post-processing: classify regulators, export tables, generate figures
#
# Input files:
#   ../4_DEG_GEO/GEO_ovarian_cancer_RRA_1row.rds   — RRA results (GEO)
#   ../1_DGE_AE/DE_full_OVARY_DESeq2_AE.rds        — DESeq2 results (TCGA vs AE)
#   ../0_DGE_GTEx/DE_full_OVARY_DESeq2_GTEx.rds    — DESeq2 results (TCGA vs GTEx)
#   ov_tcga_vst.tsv                                 — TCGA VST expression matrix
#   cancer_ovary_network_300bt_p1e-8.txt            — ARACNe-AP network
#
# Output files:
#   RRA_vs_AEDE_stat_combined.png
#   RRA_vs_GTExDE_stat_combined.png
#   core_mra.tsv
#   top_15n15_mrs.tsv
#   MRA_signature_vs_absNES_sizeFDR.png / .pdf
#   ovary_mra.RData
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(vroom)
  library(viper)
  library(parallel)
  library(readr)
})


# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

geo_rra_1row     <- readRDS("../4_DEG_GEO/GEO_ovarian_cancer_RRA_1row.rds")
DE_full_OVARY_AE <- readRDS("../1_DGE_AE/DE_full_OVARY_DESeq2_AE.rds")
DE_full_OVARY_GTEX <- readRDS("../0_DGE_GTEx/DE_full_OVARY_DESeq2_GTEx.rds")


# =============================================================================
# SECTION 1: Cross-source concordance — RRA (GEO) vs DESeq2 stat
# =============================================================================

# -----------------------------------------------------------------------------
# 1.1 Utility functions
# -----------------------------------------------------------------------------

# Convert RRA output to a signed Z-score:
#   score_rra = sign(direction) * qnorm(1 - p_adj / 2)
# Treats p_adj as a two-sided p-value; clamps to (eps, 1-eps) to avoid
# infinite Z-scores.
make_rra_signed_z <- function(tbl,
                               gene_col = "gene",
                               padj_col = "p_adj",
                               dir_col  = "direction",
                               eps      = 1e-300,
                               dedup    = TRUE) {

  stopifnot(all(c(gene_col, padj_col, dir_col) %in% colnames(tbl)))

  out <- tbl %>%
    dplyr::transmute(
      gene      = as.character(.data[[gene_col]]),
      p_adj     = pmin(pmax(as.numeric(.data[[padj_col]]), eps), 1 - eps),
      direction = as.character(.data[[dir_col]])
    ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene), direction %in% c("up", "down")) %>%
    dplyr::mutate(
      sign_dir  = ifelse(direction == "up", 1, -1),
      z_abs     = qnorm(1 - p_adj / 2),
      score_rra = sign_dir * z_abs
    ) %>%
    dplyr::select(gene, score_rra, direction, p_adj)

  if (dedup) {
    out <- out %>%
      dplyr::arrange(gene, dplyr::desc(abs(score_rra))) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
  }

  out
}

# Extract signed DESeq2 Wald statistic; deduplicate by keeping max |stat|
make_de_stat <- function(tbl,
                          gene_col = "Symbol",
                          stat_col = "stat") {

  stopifnot(all(c(gene_col, stat_col) %in% colnames(tbl)))

  tbl %>%
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      stat = as.numeric(.data[[stat_col]])
    ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene), is.finite(stat)) %>%
    dplyr::arrange(gene, dplyr::desc(abs(stat))) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
}

# Robust Spearman correlation (pairwise complete observations)
spearman_rho <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}


# -----------------------------------------------------------------------------
# 1.2 Main comparison function: RRA vs DESeq2 stat
# -----------------------------------------------------------------------------

# Computes global rank correlation and robust subset correlation between
# the RRA signed Z-score and the DESeq2 Wald statistic. Returns merged
# data and combined figure (robust scatter + global rank-rank inset).

compare_rra_vs_deseq2stat <- function(de_tbl,
                                      geo_rra_all,
                                      label = "AE-DE",
                                      top_n = 2000,
                                      eps = 1e-300) {

  # ---- Required packages ----
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
  })

  # Armar tablas limpias
  rra <- make_rra_signed_z(geo_rra_all, eps = eps)
  de  <- make_de_stat(de_tbl)

  # Merge por gen
  df <- inner_join(rra, de, by = "gene") %>%
    mutate(
      direction = factor(direction, levels = c("down", "up")),
      # NUEVO eje X (solo para visualización)
      x_logp_signed = sign(score_rra) * (-log10(p_adj))
    )

  # ----------- Global correlation (rank(|.|) vs rank(|.|)) -----------
  r1 <- rank(-abs(df$score_rra), ties.method = "average")
  r2 <- rank(-abs(df$stat),      ties.method = "average")
  rho_global <- spearman_rho(r1, r2)

  # ----------- Robust subset -----------
  df_sig <- df %>%
    filter(score_rra != 0) %>%
    arrange(desc(abs(score_rra)))

  df_top <- df_sig %>% slice_head(n = min(top_n, nrow(df_sig)))

  rho_top <- spearman_rho(df_top$score_rra, df_top$stat)

  # ----------- Etiquetas: top 10 genes extremos por lado (en df_top) -----------
  label_up <- df_top %>%
    filter(direction == "up") %>%
    arrange(desc(abs(score_rra))) %>%
    slice_head(n = 10)

  label_down <- df_top %>%
    filter(direction == "down") %>%
    arrange(desc(abs(score_rra))) %>%
    slice_head(n = 10)

  label_df <- bind_rows(label_up, label_down)

  # ----------- Figures -----------
  dir_cols <- c("up" = "red3", "down" = "dodgerblue3")

  # Panel global (inset): rank(|RRA|) vs rank(|DE stat|)
  p_global <- ggplot(df, aes(
    x = rank(-abs(score_rra), ties.method = "average"),
    y = rank(-abs(stat),      ties.method = "average"),
    color = direction
  )) +
    geom_point(alpha = 0.18, size = 0.6) +
    scale_color_manual(values = dir_cols, drop = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank(|RRA|)", y = paste0("Rank(|", label, " stat|)")) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "none",
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(2, 2, 2, 2),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 7)
    )

  # Panel robusto (principal): signed -log10(p_adj) vs stat + etiquetas
  p_robust <- ggplot(df_top, aes(x = x_logp_signed, y = stat, color = direction)) +
    geom_point(alpha = 0.35, size = 1.2) +
    scale_color_manual(values = dir_cols, drop = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.alpha = 0.6,
      show.legend = FALSE
    ) +
    labs(
      title    = "Genes with strongest RRA evidence (top |RRA|)",
      subtitle = sprintf(
        "%s: rho(top) = %.3f (N = %d) | rho(global ranks) = %.3f (N = %d)",
        label, rho_top, nrow(df_top), rho_global, nrow(df)
      ),
      x = expression(sign(RRA) %.% -log[10](p[adj])),
      y = paste0(label, " DESeq2 stat")
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )

  # Figura combinada: robusto + inset global
  p_combined <- p_robust +
    inset_element(
      p_global,
      left = 0.58, bottom = 0.07, right = 0.97, top = 0.42,
      align_to = "panel"
    )

  list(
    merged = df,
    robust_subset = df_top,
    label_points = label_df,
    rho_global_ranks = rho_global,
    rho_top_signed = rho_top,
    n_global = nrow(df),
    n_robust = nrow(df_top),
    plot_global = p_global,
    plot_robust = p_robust,
    plot_combined = p_combined
  )
}



# -----------------------------------------------------------------------------
# 1.3 Run and save concordance plots
# -----------------------------------------------------------------------------

res_ae_stat <- compare_rra_vs_deseq2stat(
  de_tbl      = DE_full_OVARY_AE,
  geo_rra_all = geo_rra_1row,
  label       = "AE-DE",
  top_n       = 2000
)

ggsave("RRA_vs_AEDE_stat_combined.png", res_ae_stat$plot_combined,
       width = 10, height = 8, dpi = 300)

res_gtex_stat <- compare_rra_vs_deseq2stat(
  de_tbl      = DE_full_OVARY_GTEX,
  geo_rra_all = geo_rra_1row,
  label       = "GTEx-DE",
  top_n       = 2000
)

ggsave("RRA_vs_GTExDE_stat_combined.png", res_gtex_stat$plot_combined,
       width = 10, height = 8, dpi = 300)

# Reported correlations
# AE:   rho(global ranks) = 0.133 (N = 14508) | rho(top robust) = 0.588 (N = 2000)
# GTEx: rho(global ranks) = 0.181 (N = 14478) | rho(top robust) = 0.658 (N = 2000)
cat("=== AE-DE ===\n")
cat("rho(global ranks):", round(res_ae_stat$rho_global_ranks, 3),
    " N =", res_ae_stat$n_global, "\n")
cat("rho(top robust):  ", round(res_ae_stat$rho_top_signed, 3),
    " N =", res_ae_stat$n_robust, "\n")

cat("=== GTEx-DE ===\n")
cat("rho(global ranks):", round(res_gtex_stat$rho_global_ranks, 3),
    " N =", res_gtex_stat$n_global, "\n")
cat("rho(top robust):  ", round(res_gtex_stat$rho_top_signed, 3),
    " N =", res_gtex_stat$n_robust, "\n")


# =============================================================================
# SECTION 2: Core consensus signature construction
# =============================================================================

# -----------------------------------------------------------------------------
# 2.1 Utility functions for signature building
# -----------------------------------------------------------------------------

# Standard Z-score normalization
zscore <- function(x) {
  x <- as.numeric(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Elbow (knee-point) detection on a sorted absolute-value vector.
# Identifies the index of maximum curvature by finding the point of greatest
# vertical distance from the diagonal baseline (normalized rank vs value).
pick_knee <- function(x) {
  x <- sort(abs(x), decreasing = TRUE)
  x <- x[is.finite(x)]
  if (length(x) < 10) return(length(x))

  i   <- seq_along(x)
  i_n <- (i - min(i)) / (max(i) - min(i))
  x_n <- (x - min(x)) / (max(x) - min(x))

  base <- 1 - i_n           # decreasing diagonal baseline
  d    <- x_n - base        # vertical distance above baseline
  which.max(d)
}

# Build a named numeric vector (msVIPER-compatible signature) from the
# core concordant gene table using z_mean as score.
make_signature_from_core <- function(core_sel,
                                      gene_col  = "gene",
                                      score_col = "z_mean") {
  stopifnot(gene_col  %in% colnames(core_sel))
  stopifnot(score_col %in% colnames(core_sel))

  sig <- core_sel[[score_col]]
  names(sig) <- core_sel[[gene_col]]
  sig <- sig[is.finite(sig) & !is.na(names(sig)) & names(sig) != ""]
  sort(sig, decreasing = TRUE)
}

# Build a named numeric vector from a DESeq2 results table.
# Optionally Z-score normalizes the Wald statistic to homogenize scale
# across sources before input to msVIPER.
make_signature_from_deseq <- function(de_tbl,
                                       score_col       = "stat",
                                       gene_col        = "Symbol_autho",
                                       transform       = c("none", "zscore"),
                                       dedup           = c("abs_max", "mean"),
                                       sort_decreasing = TRUE) {
  transform <- match.arg(transform)
  dedup     <- match.arg(dedup)

  stopifnot(all(c(gene_col, score_col) %in% colnames(de_tbl)))

  df <- de_tbl %>%
    dplyr::transmute(
      gene  = as.character(.data[[gene_col]]),
      score = as.numeric(.data[[score_col]])
    ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene), is.finite(score))

  if (dedup == "abs_max") {
    df <- df %>%
      dplyr::arrange(gene, dplyr::desc(abs(score))) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
  } else {
    df <- df %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  }

  if (transform == "zscore") {
    mu  <- mean(df$score, na.rm = TRUE)
    sdv <- stats::sd(df$score, na.rm = TRUE)
    if (is.finite(sdv) && sdv > 0) df$score <- (df$score - mu) / sdv
  }

  sig <- df$score
  names(sig) <- df$gene
  sig <- sig[is.finite(sig) & !is.na(names(sig)) & names(sig) != ""]
  if (sort_decreasing) sig <- sort(sig, decreasing = TRUE)

  sig
}


# -----------------------------------------------------------------------------
# 2.2 Build integrated score table
# -----------------------------------------------------------------------------

# Convert each source to a signed Z-score on a common gene set
rra_scores <- make_rra_signed_z(geo_rra_1row) %>%
  dplyr::transmute(gene, z_rra = score_rra)

ae_scores <- make_de_stat(DE_full_OVARY_AE,   gene_col = "Symbol", stat_col = "stat") %>%
  dplyr::rename(stat_ae = stat)

gtex_scores <- make_de_stat(DE_full_OVARY_GTEX, gene_col = "Symbol", stat_col = "stat") %>%
  dplyr::rename(stat_gtex = stat)

# Inner join to retain only genes measured in all three sources
m <- rra_scores %>%
  dplyr::inner_join(ae_scores,   by = "gene") %>%
  dplyr::inner_join(gtex_scores, by = "gene")

# Z-score normalize each source to a common scale before averaging
m <- m %>%
  dplyr::mutate(
    z_rra  = zscore(z_rra),
    z_ae   = zscore(stat_ae),
    z_gtex = zscore(stat_gtex)
  )

# Assign directional sign per source
m <- m %>%
  dplyr::mutate(
    s_rra  = sign(z_rra),
    s_ae   = sign(z_ae),
    s_gtex = sign(z_gtex)
  )

# Sign concordance summary
# TRUE  = consistent direction across all three sources
# FALSE = at least one source disagrees
table(m$s_rra == m$s_ae & m$s_rra == m$s_gtex)
# FALSE  TRUE
#  7784  6576


# -----------------------------------------------------------------------------
# 2.3 Core signature: concordant genes + elbow cutoff
# -----------------------------------------------------------------------------

# Retain only genes with concordant direction across RRA, AE, and GTEx
core_concordant <- m %>%
  dplyr::filter(s_rra == s_ae, s_rra == s_gtex) %>%
  dplyr::mutate(z_mean = (z_rra + z_ae + z_gtex) / 3) %>%
  dplyr::arrange(dplyr::desc(abs(z_mean)))

# Apply elbow criterion on |z_mean| to identify the inflection point
# that separates strongly concordant genes from the noise plateau
k <- pick_knee(core_concordant$z_mean)
message("Elbow cutoff: k = ", k, " genes retained in core signature")

core_sel <- core_concordant %>% dplyr::slice_head(n = k)
message("Unique genes in core signature: ", length(unique(core_sel$gene)))
# Expected: 6551

signature_core <- make_signature_from_core(core_sel)

# Quick checks
str(signature_core)
head(signature_core)


# =============================================================================
# SECTION 3: Regulon construction from ARACNe-AP network
# =============================================================================

# Load TCGA VST expression matrix (genes x samples)
# This is the same matrix used during ARACNe-AP network inference,
# ensuring consistent gene identifiers and feature space.
eset_tbl <- vroom::vroom("ov_tcga_vst.tsv")
eset <- eset_tbl %>%
  tibble::column_to_rownames(var = colnames(eset_tbl)[1]) %>%
  as.matrix()

# Convert ARACNe-AP network to VIPER-compatible regulon object
# Network parameters: 300 bootstraps, MI p-value threshold = 1e-8
network_file <- "cancer_ovary_network_300bt_p1e-8.txt"
regulon <- aracne2regulon(network_file, eset)


# =============================================================================
# SECTION 4: Master Regulator Analysis (msVIPER)
# =============================================================================

# -----------------------------------------------------------------------------
# 4.1 Permutation null model
# -----------------------------------------------------------------------------

# Shuffles gene labels of the signature vector, preserving the score
# distribution and signature size while disrupting regulator-target
# associations. Used to build the empirical null distribution for NES.
permute_signature_names <- function(sig) {
  sig_p <- sig
  names(sig_p) <- sample(names(sig_p), length(sig_p), replace = FALSE)
  sig_p
}

# Runs msVIPER with a signature-label permutation null model.
# Computes empirical two-sided p-values as:
#   p_emp = (1 + #{|NES_null| >= |NES_obs|}) / (N_perm + 1)
# Adjusts p-values by Benjamini-Hochberg (BH) procedure.
#
# Arguments:
#   signature  — named numeric vector (gene expression signature)
#   regulon    — VIPER regulon object (from aracne2regulon)
#   nperm      — number of permutations (default 10000)
#   minsize    — minimum regulon size (default 20 targets)
#   seed       — random seed for reproducibility
#   two_sided  — logical; if TRUE, uses |NES| for p-value computation
#   cores      — number of parallel workers
#   verbose    — print progress messages
msviper_perm_null <- function(signature,
                               regulon,
                               nperm     = 1000,
                               minsize   = 20,
                               seed      = 1,
                               two_sided = TRUE,
                               cores     = 1,
                               verbose   = TRUE) {

  stopifnot(is.numeric(signature), !is.null(names(signature)))
  set.seed(seed)

  # Observed msVIPER run
  obs      <- msviper(signature, regulon, minsize = minsize, verbose = FALSE)
  obs_nes  <- obs$es$nes
  tf_names <- names(obs_nes)

  if (verbose) {
    message("Observed msVIPER done. TFs evaluated: ", length(tf_names))
    message("Permutations: ", nperm, " | cores: ", cores, " | minsize: ", minsize)
  }

  # Single permutation: returns NES vector aligned to tf_names
  one_perm <- function(i) {
    sig_p <- permute_signature_names(signature)
    tmp   <- msviper(sig_p, regulon, minsize = minsize, verbose = FALSE)
    tmp$es$nes[tf_names]
  }

  # Run permutations (optionally in parallel)
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("signature", "regulon", "minsize", "tf_names",
                  "permute_signature_names"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, library(viper))
    null_list <- parallel::parLapply(cl, X = seq_len(nperm), fun = one_perm)
  } else {
    null_list <- lapply(seq_len(nperm), one_perm)
  }

  null_mat <- do.call(cbind, null_list)
  colnames(null_mat) <- paste0("perm", seq_len(nperm))
  rownames(null_mat) <- tf_names

  # Empirical p-values
  if (two_sided) {
    p_emp <- vapply(tf_names, function(tf) {
      (1 + sum(abs(null_mat[tf, ]) >= abs(obs_nes[tf]), na.rm = TRUE)) / (nperm + 1)
    }, numeric(1))
  } else {
    p_emp <- vapply(tf_names, function(tf) {
      (1 + sum(null_mat[tf, ] >= obs_nes[tf], na.rm = TRUE)) / (nperm + 1)
    }, numeric(1))
  }

  # BH correction
  padj_emp <- p.adjust(p_emp, method = "BH")

  out <- data.frame(
    TF       = tf_names,
    NES      = as.numeric(obs_nes[tf_names]),
    p_emp    = as.numeric(p_emp),
    padj_emp = as.numeric(padj_emp),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$padj_emp, -abs(out$NES)), ]

  list(
    observed = obs,
    null_nes = null_mat,
    table    = out,
    params   = list(nperm = nperm, minsize = minsize, seed = seed,
                    two_sided = two_sided)
  )
}


# -----------------------------------------------------------------------------
# 4.2 Run msVIPER for all three signatures
# -----------------------------------------------------------------------------

# Core consensus signature
mra_core <- msviper_perm_null(
  signature = signature_core,
  regulon   = regulon,
  nperm     = 10000,
  minsize   = 20,
  seed      = 42,
  two_sided = TRUE,
  cores     = 80,
  verbose   = TRUE
)

head(mra_core$table, 20)

# Independent validation signatures (z-score normalized DESeq2 stat)
# Z-score normalization applied to homogenize scale with signature_core
# and enable direct NES comparison across analyses.
sig_ae_z   <- make_signature_from_deseq(DE_full_OVARY_AE,   transform = "zscore")
sig_gtex_z <- make_signature_from_deseq(DE_full_OVARY_GTEX, transform = "zscore")

mra_ae_z <- msviper_perm_null(
  signature = sig_ae_z,
  regulon   = regulon,
  nperm     = 10000,
  minsize   = 20,
  seed      = 42,
  two_sided = TRUE,
  cores     = 80,
  verbose   = TRUE
)

mra_gtex_z <- msviper_perm_null(
  signature = sig_gtex_z,
  regulon   = regulon,
  nperm     = 10000,
  minsize   = 20,
  seed      = 42,
  two_sided = TRUE,
  cores     = 80,
  verbose   = TRUE
)


# =============================================================================
# SECTION 5: Cross-signature validation
# =============================================================================

# -----------------------------------------------------------------------------
# 5.1 NES rank concordance across signatures
# -----------------------------------------------------------------------------

# Interpretation:
#   high core–AE rho   → core signature dominated by AE contrast
#   high core–GTEx rho → genuine consensus signal
#   low AE–GTEx but high with core → core filters cohort-specific noise

cmp <- mra_core$table %>%
  dplyr::select(TF, NES_core = NES) %>%
  dplyr::inner_join(
    mra_ae_z$table   %>% dplyr::select(TF, NES_ae   = NES), by = "TF"
  ) %>%
  dplyr::inner_join(
    mra_gtex_z$table %>% dplyr::select(TF, NES_gtex = NES), by = "TF"
  )

cat("=== NES Spearman correlations across signatures ===\n")
cat("core vs AE:   ", round(cor(cmp$NES_core, cmp$NES_ae,   method = "spearman"), 3), "\n")
cat("core vs GTEx: ", round(cor(cmp$NES_core, cmp$NES_gtex, method = "spearman"), 3), "\n")
cat("AE vs GTEx:   ", round(cor(cmp$NES_ae,   cmp$NES_gtex, method = "spearman"), 3), "\n")


# -----------------------------------------------------------------------------
# 5.2 Overlap of significant TFs (padj_emp < 0.05)
# -----------------------------------------------------------------------------

sig_tf_core  <- mra_core$table  %>% dplyr::filter(padj_emp < 0.05) %>% dplyr::pull(TF)
sig_tf_ae    <- mra_ae_z$table  %>% dplyr::filter(padj_emp < 0.05) %>% dplyr::pull(TF)
sig_tf_gtex  <- mra_gtex_z$table %>% dplyr::filter(padj_emp < 0.05) %>% dplyr::pull(TF)

cat("=== Significant TF overlaps ===\n")
cat("AE ∩ GTEx:    ", length(intersect(sig_tf_ae,   sig_tf_gtex)), "\n")
cat("core ∩ AE:    ", length(intersect(sig_tf_core, sig_tf_ae)),   "\n")
cat("core ∩ GTEx:  ", length(intersect(sig_tf_core, sig_tf_gtex)), "\n")


# -----------------------------------------------------------------------------
# 5.3 Sign concordance (activator vs repressor direction)
# -----------------------------------------------------------------------------

# A TF positive in all three signatures = robust activated regulator
cat("=== NES sign concordance: core vs AE ===\n")
print(table(sign(cmp$NES_core), sign(cmp$NES_ae)))


# =============================================================================
# SECTION 6: Post-processing and output
# =============================================================================

# -----------------------------------------------------------------------------
# 6.1 Classify regulators from core MRA
# -----------------------------------------------------------------------------

core_mrs <- dplyr::as_tibble(mra_core$table) %>%
  dplyr::arrange(dplyr::desc(abs(NES))) %>%
  dplyr::mutate(rank_abs_NES = dplyr::row_number())

# "Is TF? == Yes" & "TF assessment == Known motif" in DatabaseExtract_v_1.01.txt from humantfs.ccbr.utoronto.ca
TFlist <- readLines("centroids/onlyTFs.txt")


core_mrs_processed <- core_mrs %>%
  dplyr::mutate(
    # Activity classification (FDR < 0.05, |NES| >= 2)
    direction = dplyr::case_when(
      padj_emp < 0.05 & NES >=  2 ~ "Activated",
      padj_emp < 0.05 & NES <= -2 ~ "Repressed",
      TRUE                        ~ "Not sig."
    ) %>% factor(levels = c("Activated", "Repressed", "Not sig.")),

    # TF expression score in the core signature
    signature_score = signature_core[TF],

    # Direction of TF expression in tumors
    expr_direction = dplyr::case_when(
      is.na(signature_score) ~ "Not in signature",
      signature_score > 0    ~ "Up in tumors",
      signature_score < 0    ~ "Down in tumors",
      TRUE                   ~ "Neutral"
    ),

    # Combined regulatory mode: activity + expression
    regulatory_mode = dplyr::case_when(
      is.na(signature_score)                                    ~ "Unknown (no signature score)",
      padj_emp < 0.05 & NES >=  2 & signature_score > 0        ~ "Activated & upregulated",
      padj_emp < 0.05 & NES >=  2 & signature_score < 0        ~ "Activated but downregulated",
      padj_emp < 0.05 & NES <= -2 & signature_score < 0        ~ "Repressed & downregulated",
      padj_emp < 0.05 & NES <= -2 & signature_score > 0        ~ "Repressed but upregulated",
      TRUE                                                      ~ "Other"
    ),
    # Whether the regulator is an annotated transcription factor
    isTF = TF %in% TFlist
  ) %>%
  dplyr::arrange(dplyr::desc(NES))

readr::write_tsv(core_mrs_processed, file = "core_mra.tsv")


# -----------------------------------------------------------------------------
# 6.2 MRA scatter plot: TF expression score vs |NES|
# -----------------------------------------------------------------------------

fdr_cut        <- 0.05
nes_cut        <- 2
label_n        <- 15   # top N per side to label

df_vol <- core_mrs_processed[core_mrs_processed$isTF== TRUE,] %>%
  dplyr::mutate(
    # TFs absent from signature are offset to the left (not zero)
    sig_x = ifelse(
      is.na(signature_score),
      min(signature_score, na.rm = TRUE) - 0.8,
      signature_score
    ),
    abs_NES = abs(NES),
    sig_fdr = padj_emp < fdr_cut,
    group4  = dplyr::case_when(
      sig_fdr & NES >=  nes_cut ~ "Activated (|NES|\u22652)",
      sig_fdr & NES <= -nes_cut ~ "Repressed (|NES|\u22652)",
      sig_fdr                   ~ "Significant (|NES|<2)",
      TRUE                      ~ "Not significant"
    ) %>% factor(levels = c(
      "Activated (|NES|\u22652)",
      "Repressed (|NES|\u22652)",
      "Significant (|NES|<2)",
      "Not significant"
    ))
  )

# Labels: top label_n activated + top label_n repressed by |NES| (FDR < 0.05)
df_labels <- dplyr::bind_rows(
  df_vol %>% dplyr::filter(sig_fdr, NES > 0) %>%
    dplyr::arrange(dplyr::desc(abs_NES)) %>% dplyr::slice_head(n = label_n),
  df_vol %>% dplyr::filter(sig_fdr, NES < 0) %>%
    dplyr::arrange(dplyr::desc(abs_NES)) %>% dplyr::slice_head(n = label_n)
)

readr::write_tsv(df_labels, file = "top_15n15_mrs.tsv")

p_mra <- ggplot(df_vol, aes(x = sig_x, y = abs_NES)) +
  geom_hline(yintercept = nes_cut, linetype = "dashed") +
  geom_vline(xintercept = 0,       linetype = "dashed") +
  geom_point(aes(color = group4, size = sig_fdr), alpha = 0.65) +
  ggrepel::geom_text_repel(
    data          = df_labels,
    aes(label = TF),
    size          = 3.0,
    box.padding   = 0.35,
    point.padding = 0.25,
    max.overlaps  = Inf
  ) +
  scale_color_manual(values = c(
    "Activated (|NES|\u22652)" = "red3",
    "Repressed (|NES|\u22652)" = "dodgerblue3",
    "Significant (|NES|<2)"    = "grey45",
    "Not significant"          = "grey80"
  )) +
  scale_size_manual(values = c(`TRUE` = 2.2, `FALSE` = 1.1), guide = "none") +
  annotate(
    "text",
    x = min(df_vol$sig_x, na.rm = TRUE),
    y = max(df_vol$abs_NES, na.rm = TRUE),
    label = "TF not in signature",
    hjust = 0, vjust = 1, size = 3
  ) +
  labs(
    title    = "Master Regulator Analysis: TF expression vs regulatory activity",
    subtitle = sprintf(
      "Color: strong (FDR<%.2g & |NES|\u2265%.1f) vs weak significant (FDR<%.2g & |NES|<%.1f). Labels: top %d per side (FDR<%.2g).",
      fdr_cut, nes_cut, fdr_cut, nes_cut, label_n, fdr_cut
    ),
    x     = "TF differential-expression score (signature_core)",
    y     = expression("|NES| (msVIPER)"),
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

print(p_mra)

ggsave("MRA_signature_vs_absNES_sizeFDR.png", p_mra, width = 10, height = 7, dpi = 300)
ggsave("MRA_signature_vs_absNES_sizeFDR.pdf", p_mra, width = 10, height = 7)

