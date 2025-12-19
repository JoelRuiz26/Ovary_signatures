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

out_dir <- file.path(base_dir,"2_Consensus_list")

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
        pvalue < P_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "up",
        pvalue < P_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "down",
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
      direction = tolower(direction)
    ) %>%
    filter(
      !is.na(gene),
      p_raw < P_CUTOFF,
      direction %in% c("up", "down")
    ) %>%
    distinct(gene, direction) %>%
    mutate(source = "GEO_RRA")
}

raw_sig <- deseq_to_signed(DE_raw, "DESeq2_raw")
ae_sig  <- deseq_to_signed(DE_ae,  "DESeq2_autoencoder")
geo_sig <- geo_to_signed(GEO)

# ===================== VENN LIST (GENES ONLY) ===================== #
venn_list <- list(
  DESeq2_raw         = unique(raw_sig$gene),
  DESeq2_autoencoder = unique(ae_sig$gene),
  GEO_RRA            = unique(geo_sig$gene)
)

# ===================== VENN PLOT ===================== #
venn_plot <- ggvenn(
  venn_list,
  fill_alpha      = 0.45,
  stroke_size     = 0.9,
  set_name_size   = 5.5,
  text_size       = 5.5,
  show_percentage = TRUE,
  digits          = 1
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

out_pdf <- file.path(out_dir, "2_0_1_venn_DESeq2_raw_DESeq2_ae_GEO.pdf")
ggsave(out_pdf, venn_plot, width = 8, height = 6.5)

# ===================== CONSENSUS GENE LIST (>=2 LISTS, SAME DIRECTION) ===================== #
consensus_genes <- bind_rows(raw_sig, ae_sig, geo_sig) %>%
  group_by(gene) %>%
  summarise(
    n_sources = n_distinct(source),
    n_up   = sum(direction == "up"),
    n_down = sum(direction == "down"),
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
