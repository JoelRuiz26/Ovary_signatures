# =========================================
# Directional meta-analysis of DE results using RobustRankAggreg
# (GTEx-normal controls)
# =========================================
# This script:
#  1) Loads gene-level DE results from edgeR, DESeq2, and limma.
#  2) Merges log2FC and p-values into a common table.
#  3) Performs rank-based meta-analysis separately for
#     up-regulated and down-regulated genes using RobustRankAggreg (RRA).
#  4) Builds a meta table (meta_all) with:
#       - average log2FC across methods
#       - direction-aware meta-analytic p-value and FDR
#       - number of individually significant methods.
#  5) Defines a consensus set of meta-significant, directionally consistent genes
#     with minimal effect size (meta_sig).
# =========================================

# Load Required Libraries
# =========================================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)
library(RobustRankAggreg)  # Rank-based meta-analysis (Kolde et al., Bioinformatics 2012)

setwd("~/Ovary_signatures/2_DGE_signature_autoEncoder/")

# =========================================
# 1. Load DE results from the three methods
# =========================================
edge   <- readRDS("2_1_Output_gtex_rds/DE_OVARY_ALL_EdgeR.rds")
deseq2 <- readRDS("2_1_Output_gtex_rds/DE_OVARY_ALL_DESeq2.rds")
limma  <- readRDS("2_1_Output_gtex_rds/DE_OVARY_ALL_limma.rds")

# =========================================
# 2. Merge DE results (logFC + p-values)
# =========================================
# For each method, we keep:
#   - log2FoldChange
#   - raw p-value
#   - adjusted p-value (FDR)

meta_tmp <- edge %>%
  dplyr::select(
    Symbol,
    log2FC_edgeR = log2FoldChange,
    p_edgeR      = pvalue,
    padj_edgeR   = padj
  ) %>%
  dplyr::inner_join(
    deseq2 %>%
      dplyr::select(
        Symbol,
        log2FC_DESeq2 = log2FoldChange,
        p_DESeq2      = pvalue,
        padj_DESeq2   = padj
      ),
    by = "Symbol"
  ) %>%
  dplyr::inner_join(
    limma %>%
      dplyr::select(
        Symbol,
        log2FC_limma = log2FoldChange,
        p_limma      = pvalue,
        padj_limma   = padj
      ),
    by = "Symbol"
  )

eps     <- .Machine$double.xmin
N_genes <- nrow(meta_tmp)

# =========================================
# 3. Directional rank-based meta-analysis with RobustRankAggreg
# =========================================
# We perform RRA separately for:
#   - Up-regulated genes (logFC > 0)
#   - Down-regulated genes (logFC < 0)
# and later assign the meta p-value according to the direction
# of the average log2FC.

# 3.1 Ranked lists for UP-regulated genes
rank_edgeR_up <- meta_tmp %>%
  dplyr::filter(log2FC_edgeR > 0) %>%
  arrange(p_edgeR) %>%
  pull(Symbol)

rank_DESeq2_up <- meta_tmp %>%
  dplyr::filter(log2FC_DESeq2 > 0) %>%
  arrange(p_DESeq2) %>%
  pull(Symbol)

rank_limma_up <- meta_tmp %>%
  dplyr::filter(log2FC_limma > 0) %>%
  arrange(p_limma) %>%
  pull(Symbol)

glist_up <- list(
  edgeR  = rank_edgeR_up,
  DESeq2 = rank_DESeq2_up,
  limma  = rank_limma_up
)

# 3.2 Ranked lists for DOWN-regulated genes
rank_edgeR_down <- meta_tmp %>%
  dplyr::filter(log2FC_edgeR < 0) %>%
  arrange(p_edgeR) %>%
  pull(Symbol)

rank_DESeq2_down <- meta_tmp %>%
  dplyr::filter(log2FC_DESeq2 < 0) %>%
  arrange(p_DESeq2) %>%
  pull(Symbol)

rank_limma_down <- meta_tmp %>%
  dplyr::filter(log2FC_limma < 0) %>%
  arrange(p_limma) %>%
  pull(Symbol)

glist_down <- list(
  edgeR  = rank_edgeR_down,
  DESeq2 = rank_DESeq2_down,
  limma  = rank_limma_down
)

# 3.3 Run RRA for UP-regulated genes
rra_up <- aggregateRanks(
  glist  = glist_up,
  N      = N_genes,
  method = "RRA",
  full   = TRUE,
  exact  = TRUE
)

colnames(rra_up) <- c("Symbol", "score_up")

rra_up <- rra_up %>%
  dplyr::mutate(
    pvalue_meta_up = pmax(score_up, eps)
  ) %>%
  dplyr::select(Symbol, pvalue_meta_up)

# 3.4 Run RRA for DOWN-regulated genes
rra_down <- aggregateRanks(
  glist  = glist_down,
  N      = N_genes,
  method = "RRA",
  full   = TRUE,
  exact  = TRUE
)

colnames(rra_down) <- c("Symbol", "score_down")

rra_down <- rra_down %>%
  dplyr::mutate(
    pvalue_meta_down = pmax(score_down, eps)
  ) %>%
  dplyr::select(Symbol, pvalue_meta_down)

# =========================================
# 4. Build meta_all table (direction-aware meta summary)
# =========================================
# For each gene:
#   - compute average log2FC across methods
#   - assign regulation (up/down)
#   - attach directional meta p-value (up or down RRA)
#   - apply a single BH correction across all genes

meta_all <- meta_tmp %>%
  dplyr::mutate(
    log2FC_mean = (log2FC_edgeR + log2FC_DESeq2 + log2FC_limma) / 3,
    regulation  = dplyr::if_else(log2FC_mean > 0, "up", "down"),
    n_sig_methods = (padj_edgeR  < 0.01) +
      (padj_DESeq2 < 0.01) +
      (padj_limma  < 0.01)
  ) %>%
  dplyr::left_join(rra_up,   by = "Symbol") %>%
  dplyr::left_join(rra_down, by = "Symbol") %>%
  dplyr::mutate(
    pvalue_meta_raw = dplyr::case_when(
      log2FC_mean > 0  ~ pvalue_meta_up,
      log2FC_mean < 0  ~ pvalue_meta_down,
      TRUE             ~ NA_real_
    ),
    pvalue_meta = pmax(dplyr::coalesce(pvalue_meta_raw, 1), eps),
    padj        = p.adjust(pvalue_meta, method = "BH"),
    padj        = pmax(padj, eps)
  ) %>%
  dplyr::select(
    Symbol,
    log2FoldChange = log2FC_mean,
    regulation,
    pvalue_meta,
    padj,
    n_sig_methods,
    log2FC_edgeR,  p_edgeR,  padj_edgeR,
    log2FC_DESeq2, p_DESeq2, padj_DESeq2,
    log2FC_limma,  p_limma,  padj_limma
  )

# Save full meta-analysis table
saveRDS(
  meta_all,
  "2_1_Output_gtex_rds/2_2_Meta_DE_OVARY_ALL_RobustRankAggreg_directional.rds"
)

# Optional comparison table
meta_all_compare <- meta_all %>%
  dplyr::select(
    Symbol,
    pvalue_meta,
    padj,
    n_sig_methods,
    p_edgeR,
    p_DESeq2,
    p_limma
  )

# =========================================
# 5. Define meta-significant & consistent genes
# =========================================
# Criteria:
#   - Same sign of log2FC across edgeR, DESeq2, and limma
#   - Direction-aware meta FDR (padj) <= 0.01
#   - At least 2 of the 3 methods individually significant (FDR < 0.01)
#   - Minimum effect size: |log2FoldChange| >= 1

meta_sig <- meta_all %>%
  dplyr::filter(
    sign(log2FC_edgeR)   == sign(log2FC_DESeq2),
    sign(log2FC_edgeR)   == sign(log2FC_limma),
    padj <= 0.01,
    n_sig_methods >= 2,
    abs(log2FoldChange) >= 1
  )

saveRDS(
  meta_sig,
  "2_1_Output_gtex_rds/2_2_Meta_DE_OVARY_significant_RobustRankAggreg_directional.rds"
)
