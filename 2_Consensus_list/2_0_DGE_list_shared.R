#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggvenn)
  library(ggplot2)
  library(readr)
})

options(stringsAsFactors = FALSE)

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
# Circle membership: p < 0.01 (ALL significant genes)
P_CUTOFF   <- 0.01
# Direction is defined only if |log2FC| > 1 for DESeq2
LFC_CUTOFF <- 1

# ===================== HELPERS ===================== #
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# Return ALL significant genes (p < cutoff), gene-only (for circle size)
deseq_sig_genes_only <- function(df) {
  df %>%
    mutate(gene = clean_gene(Symbol)) %>%
    filter(!is.na(gene), !is.na(pvalue), pvalue < P_CUTOFF) %>%
    distinct(gene) %>%
    pull(gene)
}

geo_sig_genes_only <- function(df) {
  df %>%
    mutate(gene = clean_gene(Gene.symbol)) %>%
    filter(!is.na(gene), !is.na(p_raw), p_raw < P_CUTOFF) %>%
    distinct(gene) %>%
    pull(gene)
}

# Return significant genes WITH a well-defined direction (used for concordance logic)
deseq_to_signed <- function(df, source_name) {
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
    distinct(gene, direction) %>%
    mutate(source = source_name)
}

geo_to_signed <- function(df) {
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

# Build:
# (A) full significant sets (circle membership)
raw_all <- deseq_sig_genes_only(DE_raw)
ae_all  <- deseq_sig_genes_only(DE_ae)
geo_all <- geo_sig_genes_only(GEO)

# (B) directed significant tables (for concordance)
raw_sig <- deseq_to_signed(DE_raw, "TCGA (ovary_GTEX_control)")
ae_sig  <- deseq_to_signed(DE_ae,  "TCGA (Reference_control)")
geo_sig <- geo_to_signed(GEO)

# ===================== CONSENSUS GENE LIST (WRITE OUTPUT) ===================== #
# Concordant direction in >=2 sources (up in >=2 OR down in >=2)
consensus_genes <- bind_rows(raw_sig, ae_sig, geo_sig) %>%
  distinct(gene, source, direction) %>%   # one vote per source
  group_by(gene) %>%
  summarise(
    n_sources = n_distinct(source),
    n_up      = n_distinct(source[direction == "up"]),
    n_down    = n_distinct(source[direction == "down"]),
    consensus_direction = case_when(
      n_up   >= 2 ~ "up",
      n_down >= 2 ~ "down",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(n_sources >= 2, !is.na(consensus_direction)) %>%
  arrange(desc(n_sources), gene)

out_genes <- file.path(out_dir, "2_0_1_genes_concordant.tsv")
write_tsv(consensus_genes, out_genes)

# ===================== VENN SETS: FULL CIRCLES, CONCORDANT OVERLAPS ===================== #
# Trick for ggvenn:
# - Keep concordant genes as the SAME ID across sources => they intersect.
# - Rename non-concordant (or no-direction) genes per source => they stay in the circle but do NOT intersect.

# Consensus direction lookup: gene -> consensus_direction
cons_dir <- setNames(consensus_genes$consensus_direction, consensus_genes$gene)

# Per-source direction lookup (only when direction is uniquely defined)
make_dir_map <- function(sig_df) {
  sig_df %>%
    group_by(gene) %>%
    summarise(dir = ifelse(n_distinct(direction) == 1, first(direction), NA_character_), .groups = "drop") %>%
    { setNames(.$dir, .$gene) }
}

raw_dir <- make_dir_map(raw_sig)
ae_dir  <- make_dir_map(ae_sig)
geo_dir <- make_dir_map(geo_sig)

# Build plot IDs:
# If gene is concordant AND this source matches the consensus direction => keep "GENE"
# Else => rename to "GENE__<TAG>" (still counted in the circle, but won't overlap)
make_plot_ids <- function(gene_vec, dir_map, tag) {
  vapply(gene_vec, function(g) {
    cd <- unname(cons_dir[g])
    sd <- unname(dir_map[g])
    if (!is.na(cd) && !is.na(sd) && sd == cd) g else paste0(g, "__", tag)
  }, character(1))
}


raw_plot <- unique(make_plot_ids(raw_all, raw_dir, "RAW"))
ae_plot  <- unique(make_plot_ids(ae_all,  ae_dir,  "AE"))
geo_plot <- unique(make_plot_ids(geo_all, geo_dir, "GEO"))

venn_list <- list(
  "TCGA\n(ovary_GTEX_control)" = raw_plot,
  "TCGA\n(Reference_control)"  = ae_plot,
  "GEO\n(6 datasets)"          = geo_plot
)

# ===================== VENN PLOT (UNCHANGED TITLES / OUTPUT NAMES) ===================== #
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

out_pdf <- file.path(out_dir, "2_0_1_venn_TCGA_GEO.pdf")
ggsave(out_pdf, venn_plot, width = 8, height = 6.5)

cat("Done.\n")
cat("Venn diagram saved to: ", out_pdf, "\n", sep = "")
cat("Consensus gene list saved to: ", out_genes, "\n", sep = "")
