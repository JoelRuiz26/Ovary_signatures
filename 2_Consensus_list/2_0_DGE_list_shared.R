#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggvenn)
  library(ggplot2)
  library(readr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Goal
# 1) Venn diagram of significant genes:
#    - DESeq2: p < 0.01 AND |log2FC| > 1
#    - GEO RRA: p_raw < 0.01 (direction already provided)
#    (Same labels, title, colors, and output filenames.)
#
# 2) Table of genes present in >=2 sources with:
#    gene, n_sources, n_up, n_down, consensus_direction
#    consensus_direction = "up" if n_up > n_down
#                        = "down" if n_down > n_up
#                        = "none" if tie (common when gene appears in exactly 2 sources with opposite directions)
# ============================================================

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw_deseq2 = file.path(base_dir, "0_DGE_raw",         "DE_full_OVARY_DESeq2.rds"),
  ae_deseq2  = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds"),
  geo_rra    = "/STORAGE/csbig/hachepunto/adenosina/Ovary_signatures/4_DEG_GEO/GEO_ovarian_cancer_RRA_results.rds"
)

out_dir <- file.path(base_dir, "2_Consensus_list")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("File not found: ", p)
  readRDS(p)
}

DE_raw <- read_rds_safe(paths$raw_deseq2)
DE_ae  <- read_rds_safe(paths$ae_deseq2)
GEO    <- read_rds_safe(paths$geo_rra)

# ===================== THRESHOLDS ===================== #
P_CUTOFF   <- 0.01
LFC_CUTOFF <- 1

# ===================== HELPERS ===================== #
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# ===================== BUILD "SIGNED" LISTS (WHAT DEFINES SIGNIFICANCE) ===================== #
# These are the ONLY genes that will go into the Venn.
# (So Venn counts will match the logic used in the table.)

deseq_sig_signed <- function(df, source_name) {
  df %>%
    mutate(
      gene = clean_gene(Symbol),
      direction = case_when(
        !is.na(pvalue) & pvalue < P_CUTOFF & !is.na(log2FoldChange) & log2FoldChange >  LFC_CUTOFF ~ "up",
        !is.na(pvalue) & pvalue < P_CUTOFF & !is.na(log2FoldChange) & log2FoldChange < -LFC_CUTOFF ~ "down",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(gene), !is.na(direction)) %>%
    distinct(gene, direction) %>%  # keep both if a gene truly appears with both directions (rare)
    mutate(source = source_name)
}

geo_sig_signed <- function(df) {
  df %>%
    mutate(
      gene = clean_gene(Gene.symbol),
      direction = tolower(trimws(direction))
    ) %>%
    filter(
      !is.na(gene),
      !is.na(p_raw), p_raw < P_CUTOFF,
      direction %in% c("up", "down")
    ) %>%
    distinct(gene, direction) %>%
    mutate(source = "GEO")
}

raw_sig <- deseq_sig_signed(DE_raw, "TCGA (ovary_GTEX_control)")
ae_sig  <- deseq_sig_signed(DE_ae,  "TCGA (Reference_control)")
geo_sig <- geo_sig_signed(GEO)

# ===================== TABLE: GENES PRESENT IN >=2 SOURCES ===================== #
# One vote per source per gene:
# - If a gene has multiple rows within a source (shouldn't happen often), we collapse to a single direction vote per source.
#   If it conflicts within the same source, we set that source's vote to NA (excluded from up/down voting).
collapse_to_one_vote_per_source <- function(df) {
  df %>%
    group_by(gene, source) %>%
    summarise(
      dir = if (n_distinct(direction) == 1) first(direction) else NA_character_,
      .groups = "drop"
    ) %>%
    filter(!is.na(dir)) %>%
    rename(direction = dir)
}

votes <- bind_rows(raw_sig, ae_sig, geo_sig) %>%
  collapse_to_one_vote_per_source()

consensus_genes <- votes %>%
  group_by(gene) %>%
  summarise(
    n_sources = n_distinct(source),
    n_up      = sum(direction == "up"),
    n_down    = sum(direction == "down"),
    consensus_direction = case_when(
      n_up   > n_down ~ "up",
      n_down > n_up   ~ "down",
      TRUE            ~ "none"   # tie (e.g., 1 up vs 1 down)
    ),
    .groups = "drop"
  ) %>%
  filter(n_sources >= 2) %>%
  arrange(desc(n_sources), gene)

# Output table (KEEP SAME FILENAME)
out_genes <- file.path(out_dir, "2_0_1_genes_concordant.tsv")
write_tsv(consensus_genes, out_genes)

# ===================== VENN: SIGNIFICANT GENES ONLY (MATCHING THE SAME RULES) ===================== #
# Venn uses only presence/absence of gene symbols in each significant list.
# (This will now be consistent with the "votes" universe above.)

venn_list <- list(
  "TCGA\n(ovary_GTEX_control)" = sort(unique(raw_sig$gene)),
  "TCGA\n(Reference_control)"  = sort(unique(ae_sig$gene)),
  "GEO\n(6 datasets)"          = sort(unique(geo_sig$gene))
)

venn_plot <- ggvenn(
  venn_list,
  fill_color      = c("#7A8DB8", "#F1E27C", "#7FC97F"),
  fill_alpha      = 0.55,
  stroke_size     = 1,
  set_name_size   = 5,
  text_size       = 6,
  show_percentage = FALSE
) +
  labs(
    title = "Ovarian cancer: significant genes across databases"
  ) +
  coord_fixed() +
  theme_classic(base_size = 14) +
  theme(
    axis.line   = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    plot.title  = element_text(face = "bold", hjust = 0.5, margin = margin(b = 18)),
    plot.margin = margin(15, 20, 15, 20)
  )

# Output Venn (KEEP SAME FILENAME)
out_pdf <- file.path(out_dir, "2_0_1_venn_TCGA_GEO.pdf")
ggsave(out_pdf, venn_plot, width = 8, height = 6.5)

cat("Done.\n")
cat("Venn diagram saved to: ", out_pdf, "\n", sep = "")
cat("Consensus gene list saved to: ", out_genes, "\n", sep = "")
