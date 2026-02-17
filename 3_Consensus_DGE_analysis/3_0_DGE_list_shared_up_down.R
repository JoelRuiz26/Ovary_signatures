#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggvenn)
  library(ggplot2)
  library(readr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Outputs
# 1) Two Venn diagrams of significant genes:
#    - UP genes:   TCGA1 vs TCGA2 vs GEO
#    - DOWN genes: TCGA1 vs TCGA2 vs GEO
#
# 2) Table of consensus genes defined BY DIRECTION (no voting):
#    gene, direction, n_sources, sources
#    where n_sources >= 2 within that direction
# ============================================================

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw_deseq2 = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae_deseq2  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  geo_rra    = file.path(base_dir, "2_DEG_GEO",  "GEO_ovarian_cancer_RRA_results.rds")
)

out_dir <- file.path(base_dir, "3_Consensus_DGE_analysis")
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

# ===================== BUILD "SIGNED" LISTS (SIGNIFICANT UP/DOWN) ===================== #
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
    distinct(gene, direction) %>%
    mutate(source = source_name)
}

geo_sig_signed <- function(df) {
  df %>%
    mutate(
      gene = clean_gene(gene),
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

# ===================== COLLAPSE TO ONE VOTE PER SOURCE PER GENE+DIR ===================== #
# (Protección contra duplicados raros dentro de una fuente)
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

# ===================== CONSENSUS (NO VOTING UP vs DOWN) ===================== #
# Consenso UP: genes "up" presentes en >=2 fuentes
# Consenso DOWN: genes "down" presentes en >=2 fuentes
consensus_by_dir <- votes %>%
  group_by(gene, direction) %>%
  summarise(
    n_sources = n_distinct(source),
    sources   = paste(sort(unique(source)), collapse = "; "),
    .groups   = "drop"
  ) %>%
  filter(n_sources >= 2) %>%
  arrange(direction, desc(n_sources), gene)

# Output table (KEEP SAME FILENAME)
out_genes <- file.path(out_dir, "3_0_1_genes_concordant.tsv")
write_tsv(consensus_by_dir, out_genes)

# ===================== VENNs (TWO: UP and DOWN) ===================== #
make_venn <- function(venn_list, title_text) {
  ggvenn(
    venn_list,
    fill_color      = c("#7A8DB8", "#F1E27C", "#7FC97F"),
    fill_alpha      = 0.55,
    stroke_size     = 1,
    set_name_size   = 5,
    text_size       = 6,
    show_percentage = FALSE
  ) +
    labs(title = title_text) +
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
}

# UP lists
venn_list_up <- list(
  "TCGA\n(ovary_GTEX_control)" = sort(unique(raw_sig$gene[raw_sig$direction == "up"])),
  "TCGA\n(Reference_control)"  = sort(unique(ae_sig$gene[ae_sig$direction == "up"])),
  "GEO\n(6 datasets)"          = sort(unique(geo_sig$gene[geo_sig$direction == "up"]))
)

venn_up <- make_venn(
  venn_list_up,
  "Overexpressed"
)

out_pdf_up <- file.path(out_dir, "3_0_1_venn_TCGA_GEO_UP.pdf")
out_png_up <- file.path(out_dir, "3_0_1_venn_TCGA_GEO_UP.png")
ggsave(out_pdf_up, venn_up, width = 8, height = 6.5)
ggsave(out_png_up, venn_up, width = 8, height = 6.5, dpi = 600)

# DOWN lists
venn_list_down <- list(
  "TCGA\n(ovary_GTEX_control)" = sort(unique(raw_sig$gene[raw_sig$direction == "down"])),
  "TCGA\n(Reference_control)"  = sort(unique(ae_sig$gene[ae_sig$direction == "down"])),
  "GEO\n(6 datasets)"          = sort(unique(geo_sig$gene[geo_sig$direction == "down"]))
)

venn_down <- make_venn(
  venn_list_down,
  "Underexpressed genes"
)

out_pdf_down <- file.path(out_dir, "3_0_1_venn_TCGA_GEO_DOWN.pdf")
out_png_down <- file.path(out_dir, "3_0_1_venn_TCGA_GEO_DOWN.png")
ggsave(out_pdf_down, venn_down, width = 8, height = 6.5)
ggsave(out_png_down, venn_down, width = 8, height = 6.5, dpi = 600)

cat("Done.\n")
cat("UP Venn saved to:   ", out_pdf_up, "\n", sep = "")
cat("DOWN Venn saved to: ", out_pdf_down, "\n", sep = "")
cat("Consensus-by-direction table saved to: ", out_genes, "\n", sep = "")
