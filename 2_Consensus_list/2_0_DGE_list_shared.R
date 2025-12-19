#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggvenn)
  library(ggplot2)
  library(stringr)
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
# Per your requirement: significance is based on *p-value* (not padj)
P_CUTOFF   <- 0.01
# Direction for DESeq2 is defined only when |log2FC| > 1
LFC_CUTOFF <- 1

# ===================== HELPERS ===================== #
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# ===================== PROCESS EACH SOURCE ===================== #
deseq_to_signed <- function(df, source_name) {
  if (!("Symbol" %in% colnames(df))) stop("DESeq2 table missing column: Symbol")
  if (!("pvalue" %in% colnames(df))) stop("DESeq2 table missing column: pvalue")
  if (!("log2FoldChange" %in% colnames(df))) stop("DESeq2 table missing column: log2FoldChange")
  
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
  if (!("Gene.symbol" %in% colnames(df))) stop("GEO table missing column: Gene.symbol")
  if (!("p_raw" %in% colnames(df))) stop("GEO table missing column: p_raw")
  if (!("direction" %in% colnames(df))) stop("GEO table missing column: direction")
  
  df %>%
    mutate(
      gene = clean_gene(Gene.symbol),
      direction = tolower(trimws(direction)),
      direction = ifelse(direction %in% c("up", "down"), direction, NA_character_)
    ) %>%
    filter(
      !is.na(gene),
      !is.na(p_raw), p_raw < P_CUTOFF,
      !is.na(direction)
    ) %>%
    distinct(gene, direction) %>%
    mutate(source = "GEO_RRA")
}

# IMPORTANT: use the set names you want to show in the Venn figure
raw_sig <- deseq_to_signed(DE_raw, "ovary_gtex_control")
ae_sig  <- deseq_to_signed(DE_ae,  "Reference_control")
geo_sig <- geo_to_signed(GEO)

# ===================== VENN LIST (GENES ONLY) ===================== #
# Venn is presence/absence only (ignores direction)
venn_list <- list(
  ovary_gtex_control = unique(raw_sig$gene),
  Reference_control  = unique(ae_sig$gene),
  GEO_RRA            = unique(geo_sig$gene)
)

# ===================== VENN PLOT (NO PERCENTAGES) ===================== #
venn_plot <- ggvenn(
  venn_list,
  fill_alpha      = 0.45,
  stroke_size     = 0.9,
  set_name_size   = 5.5,
  text_size       = 5.5,
  show_percentage = FALSE  # <<< only counts
) +
  labs(
    title = "Ovarian cancer: significant genes across datasets",
    subtitle = paste0(
      "p < ", P_CUTOFF,
      " | DESeq2 direction: |log2FC| > ", LFC_CUTOFF,
      " | GEO direction from RRA"
    )
  ) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  theme(
    axis.line  = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

out_pdf <- file.path(out_dir, "2_0_1_venn_ovary_gtex_reference_GEO.pdf")
ggsave(out_pdf, venn_plot, width = 8, height = 6.5)

# ===================== CONSENSUS GENE LIST (>=2 LISTS, SAME DIRECTION) ===================== #
# Rule:
# - gene must be significant (p < 0.01) with defined direction in >= 2 sources
# - concordant direction in >= 2 sources: (n_up >= 2) OR (n_down >= 2)
consensus_genes <- bind_rows(raw_sig, ae_sig, geo_sig) %>%
  group_by(gene) %>%
  summarise(
    n_sources = n_distinct(source),
    n_up      = sum(direction == "up"),
    n_down    = sum(direction == "down"),
    consensus_direction = case_when(
      n_up   >= 2 ~ "up",
      n_down >= 2 ~ "down",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_sources >= 2,
    !is.na(consensus_direction)
  ) %>%
  arrange(desc(n_sources), gene)

out_genes <- file.path(out_dir, "2_0_1_genes_concordant.tsv")
write_tsv(consensus_genes, out_genes)

cat("Done.\n")
cat("Venn diagram saved to: ", out_pdf, "\n", sep = "")
cat("Consensus gene list saved to: ", out_genes, "\n", sep = "")
