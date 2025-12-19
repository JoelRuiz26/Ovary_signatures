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
# Significance: p < 0.01 (as requested)
P_CUTOFF   <- 0.01
# Direction is defined only when |log2FC| > 1 for DESeq2
LFC_CUTOFF <- 1

# ===================== HELPERS ===================== #
clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# ===================== BUILD PER-SOURCE SIGNIFICANT LISTS ===================== #
# IMPORTANT: we build TWO representations per source:
# (A) gene-only significant set (p < 0.01) -> used for the outer circles (all significant genes)
# (B) gene+direction significant set (p < 0.01 AND |log2FC|>1 / GEO direction) -> used to compute concordance

deseq_sig_genes_only <- function(df) {
  df %>%
    mutate(gene = clean_gene(Symbol)) %>%
    filter(!is.na(gene), !is.na(pvalue), pvalue < P_CUTOFF) %>%
    distinct(gene) %>%
    pull(gene)
}

deseq_sig_with_direction <- function(df, source_name) {
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

geo_sig_genes_only <- function(df) {
  df %>%
    mutate(gene = clean_gene(Gene.symbol)) %>%
    filter(!is.na(gene), !is.na(p_raw), p_raw < P_CUTOFF) %>%
    distinct(gene) %>%
    pull(gene)
}

geo_sig_with_direction <- function(df) {
  df %>%
    mutate(
      gene = clean_gene(Gene.symbol),
      direction = tolower(trimws(direction))
    ) %>%
    filter(!is.na(gene), !is.na(p_raw), p_raw < P_CUTOFF, direction %in% c("up", "down")) %>%
    distinct(gene, direction) %>%
    mutate(source = "GEO")
}

# --- (A) gene-only significant sets (these define the FULL circles) ---
raw_genes_all <- deseq_sig_genes_only(DE_raw)
ae_genes_all  <- deseq_sig_genes_only(DE_ae)
geo_genes_all <- geo_sig_genes_only(GEO)

# --- (B) gene+direction significant tables (these define CONCORDANCE) ---
raw_signed <- deseq_sig_with_direction(DE_raw, "TCGA (ovary_GTEX_control)")
ae_signed  <- deseq_sig_with_direction(DE_ae,  "TCGA (Reference_control)")
geo_signed <- geo_sig_with_direction(GEO)

# ===================== CONSENSUS TABLE (CONCORDANT DIRECTION IN >=2 SOURCES) ===================== #
# One vote per source (prevents duplicated symbols from inflating up/down counts)
consensus_genes <- bind_rows(raw_signed, ae_signed, geo_signed) %>%
  distinct(gene, source, direction) %>%
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

# Output 1 (keep the SAME filename)
out_genes <- file.path(out_dir, "2_0_1_genes_concordant.tsv")
write_tsv(consensus_genes, out_genes)

# ===================== VENN (FULL SETS, BUT INTERSECTIONS = CONCORDANT) ===================== #
# Requirement:
# - Each circle should include ALL significant genes (p < 0.01) for that source.
# - Overlaps should count only genes that are concordant in direction (>=2 sources).
#
# A standard Venn cannot apply different rules to "only-in-set" vs "overlap" regions automatically.
# So we enforce the overlap rule by:
# - Keeping the circles as ALL significant genes
# - Removing "discordant-overlap genes" from the intersections by assigning those genes to ONLY ONE set.
#
# Practically: for any gene that appears in >=2 circles BUT is NOT in the concordant set,
# we assign it to a single circle (priority: raw -> ae -> geo), so it does not contribute to overlaps.

cons_set <- unique(consensus_genes$gene)

raw_set <- unique(raw_genes_all)
ae_set  <- unique(ae_genes_all)
geo_set <- unique(geo_genes_all)

# Genes that would create overlaps but are NOT concordant
overlap_any <- Reduce(intersect, list(raw_set, ae_set)) %>% union(Reduce(intersect, list(raw_set, geo_set))) %>% union(Reduce(intersect, list(ae_set, geo_set)))
discordant_overlap <- setdiff(overlap_any, cons_set)

# Remove discordant-overlap genes from AE and GEO if they are present in RAW (so they become RAW-only)
ae_set  <- setdiff(ae_set,  intersect(discordant_overlap, raw_set))
geo_set <- setdiff(geo_set, intersect(discordant_overlap, raw_set))

# Remove remaining discordant-overlap genes from GEO if they are present in AE (so they become AE-only)
geo_set <- setdiff(geo_set, intersect(discordant_overlap, ae_set))

# Now the ONLY genes that can remain in overlaps are concordant ones (cons_set)
# (Concordant genes remain in all sets where they are significant.)

venn_list <- list(
  "TCGA\n(ovary_GTEX_control)" = sort(raw_set),
  "TCGA\n(Reference_control)"  = sort(ae_set),
  "GEO\n(6 datasets)"          = sort(geo_set)
)

# ===================== VENN PLOT ===================== #
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

# Output 2 (keep the SAME filename)
out_pdf <- file.path(out_dir, "2_0_1_venn_TCGA_GEO.pdf")
ggsave(out_pdf, venn_plot, width = 8, height = 6.5)

cat("Done.\n")
cat("Venn diagram saved to: ", out_pdf, "\n", sep = "")
cat("Consensus gene list saved to: ", out_genes, "\n", sep = "")
