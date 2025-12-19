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
P_CUTOFF   <- 0.01
LFC_CUTOFF <- 1

# ===================== HELPERS ===================== #
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# ===================== PROCESS EACH SOURCE ===================== #
deseq_to_signed <- function(df, source_name) {
  df %>%
    mutate(
      gene = clean_gene(Symbol),
      direction = case_when(
        !is.na(pvalue) & pvalue < P_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "up",
        !is.na(pvalue) & pvalue < P_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "down",
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
      p_raw < P_CUTOFF,
      direction %in% c("up", "down")
    ) %>%
    distinct(gene, direction) %>%
    mutate(source = "GEO")
}

raw_sig <- deseq_to_signed(DE_raw, "TCGA (ovary_GTEX_control)")
ae_sig  <- deseq_to_signed(DE_ae,  "TCGA (Reference_control)")
geo_sig <- geo_to_signed(GEO)

# ===================== VENN LIST ===================== #
venn_list <- list(
  "TCGA\n(ovary_GTEX_control)" = unique(raw_sig$gene),
  "TCGA\n(Reference_control)"  = unique(ae_sig$gene),
  "GEO"                        = unique(geo_sig$gene)
)


# ===================== VENN PLOT ===================== #
venn_plot <- ggvenn(
  venn_list,
  fill_color     = c("#7A8DB8", "#F1E27C", "#7FC97F"),  # paper-grade palette
  fill_alpha     = 0.55,
  stroke_size    = 1,
  set_name_size  = 6,
  text_size      = 6,
  show_percentage = FALSE
) +
  labs(
    title = "Ovarian cancer: significant genes across datasets"
  ) +
  coord_fixed() +
  theme_classic(base_size = 14) +
  theme(
    axis.line  = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(15, 20, 15, 20)
  )

out_pdf <- file.path(out_dir, "2_0_1_venn_TCGA_GEO.pdf")
ggsave(out_pdf, venn_plot, width = 8, height = 6.5)

# ===================== CONSENSUS GENE LIST ===================== #
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
  filter(n_sources >= 2, !is.na(consensus_direction)) %>%
  arrange(desc(n_sources), gene)

out_genes <- file.path(out_dir, "2_0_1_genes_concordant.tsv")
write_tsv(consensus_genes, out_genes)

cat("Done.\n")
cat("Venn diagram saved to: ", out_pdf, "\n", sep = "")
cat("Consensus gene list saved to: ", out_genes, "\n", sep = "")
