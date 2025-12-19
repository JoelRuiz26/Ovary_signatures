#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(readr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Goal
# - Filter DepMap CRISPR dependency table to OVARY cell lines
# - Keep only genes in your DGE_list_shared with consensus_direction != "none"
# - Add consensus_direction (up/down) as a column to the DepMap table
# - Save ONE output (no summaries, no TPM)
# ============================================================

# -------------------- Paths -------------------- #
DGE_PATH <- "/home/jruiz/Ovary_signatures/2_Consensus_list/2_0_1_genes_concordant.tsv"
OUT_DIR  <- "/home/jruiz/Ovary_signatures/4_depmap_ovary"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------- Load consensus table -------------------- #
DGE_list_shared <- vroom(DGE_PATH, show_col_types = FALSE)

# Keep only genes with defined direction (exclude ties)
cons_genes <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(!is.na(consensus_direction), consensus_direction != "none") %>%
  distinct(gene_u, consensus_direction)

cat("Consensus genes kept (direction != none): ", nrow(cons_genes), "\n", sep = "")

# -------------------- DepMap CRISPR (OVARY only) -------------------- #
crispr_ovary <- depmap_crispr() %>%
  filter(grepl("_OVARY$", cell_line)) %>%
  mutate(gene_u = toupper(gene_name)) %>%
  inner_join(cons_genes, by = "gene_u") %>%
  transmute(
    depmap_id,
    cell_line,
    gene = gene_name,
    dependency,
    consensus_direction
  )

#> length(unique(crispr_ovary$gene))
#[1] 3078

# -------------------- Save output -------------------- #
out_crispr <- file.path(OUT_DIR, "4_1_depmap_CRISPR_OVARY_genes_in_DGE_list_shared.tsv")
write_tsv(crispr_ovary, out_crispr)

cat("Saved: ", out_crispr, "\n", sep = "")
cat("Done.\n")








# ===================== PLOTS (replace your previous plot block) ===================== #
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# Remove "_OVARY" from cell line labels (for plotting)
crispr_plot <- crispr_ovary %>%
  mutate(cell = sub("_OVARY$", "", cell_line))

# Genes to highlight (vertical markers) + keep order
highlight_genes <- c("ADA", "ADORA2A", "ADORA3")

# Common limits (same scale for UP and DOWN)
lims <- c(-3.3, 2.0)

# One shared red↔blue scale (NEGATIVE = blue = more essential; POSITIVE = red)
ceres_fill_rb <- function() {
  scale_fill_gradient2(
    low = "#2C7BB6",  # blue (more negative)
    mid = "white",
    high = "#D7191C", # red (more positive)
    midpoint = 0,
    limits = lims,
    oob = scales::squish,
    name = "Dependency\n(CERES)"
  )
}

# Rank genes by mean dependency (more negative = more essential)
gene_importance <- crispr_plot %>%
  group_by(gene, consensus_direction) %>%
  summarise(
    n_cells  = n_distinct(cell),
    mean_dep = mean(dependency, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_dep)

# Rank cell lines by mean dependency across genes in a given dataframe
rank_cells_by_mean <- function(df) {
  df %>%
    group_by(cell) %>%
    summarise(cell_impact = mean(dependency, na.rm = TRUE), .groups = "drop") %>%
    arrange(cell_impact) %>%   # most negative first
    pull(cell)
}

# Heatmap settings
top_n_genes <- 150
top_n_cells <- 60   # set Inf to keep all

# Helper to compute x positions for vertical markers
vline_df <- function(levels_gene) {
  tibble(gene = intersect(highlight_genes, levels_gene)) %>%
    mutate(xi = as.numeric(factor(gene, levels = levels_gene)))
}

# --------------------- Heatmap 1: UP consensus genes --------------------- #
genes_up <- gene_importance %>%
  filter(consensus_direction == "up") %>%
  slice_head(n = top_n_genes) %>%
  pull(gene)

df_up <- crispr_plot %>%
  filter(consensus_direction == "up", gene %in% genes_up) %>%
  mutate(gene = factor(gene, levels = genes_up))

cells_up <- rank_cells_by_mean(df_up)
if (is.finite(top_n_cells)) cells_up <- head(cells_up, top_n_cells)

df_up <- df_up %>%
  filter(cell %in% cells_up) %>%
  mutate(cell = factor(cell, levels = cells_up))

p_up <- ggplot(df_up, aes(x = gene, y = cell, fill = dependency)) +
  geom_tile() +
  geom_vline(
    data = vline_df(levels(df_up$gene)),
    aes(xintercept = xi),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  ceres_fill_rb() +
  labs(
    title = "DepMap CRISPR dependency in OVARY cell lines (Up-consensus genes)",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "4_2_heatmap_depmap_OVARY_UP.pdf"),
       p_up, width = 10, height = 7.5)

# --------------------- Heatmap 2: DOWN consensus genes --------------------- #
genes_down <- gene_importance %>%
  filter(consensus_direction == "down") %>%
  slice_head(n = top_n_genes) %>%
  pull(gene)

df_down <- crispr_plot %>%
  filter(consensus_direction == "down", gene %in% genes_down) %>%
  mutate(gene = factor(gene, levels = genes_down))

cells_down <- rank_cells_by_mean(df_down)
if (is.finite(top_n_cells)) cells_down <- head(cells_down, top_n_cells)

df_down <- df_down %>%
  filter(cell %in% cells_down) %>%
  mutate(cell = factor(cell, levels = cells_down))

p_down <- ggplot(df_down, aes(x = gene, y = cell, fill = dependency)) +
  geom_tile() +
  geom_vline(
    data = vline_df(levels(df_down$gene)),
    aes(xintercept = xi),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  ceres_fill_rb() +
  labs(
    title = "DepMap CRISPR dependency in OVARY cell lines (Down-consensus genes)",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "4_3_heatmap_depmap_OVARY_DOWN.pdf"),
       p_down, width = 10, height = 7.5)

# --------------------- Key genes panel: ADA / ADORA2A / ADORA3 --------------------- #
# Order OVARY lines *per gene* by dependency (most negative = most important)
df_markers <- crispr_plot %>%
  filter(gene %in% highlight_genes) %>%
  mutate(gene = factor(gene, levels = highlight_genes)) %>%
  group_by(gene) %>%
  mutate(cell_rank = rank(dependency, ties.method = "first")) %>%
  ungroup()

# Keep x labels only once (bottom), and don’t repeat y labels in every facet
p_markers <- ggplot(df_markers, aes(x = gene, y = reorder(cell, dependency), fill = dependency)) +
  geom_tile() +
  facet_grid(. ~ gene, scales = "free_y", space = "free_y") +
  # Different palette for the callout plot (still red↔blue, but higher contrast)
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0,
    limits = lims,
    oob = scales::squish,
    name = "Dependency\n(CERES)"
  ) +
  labs(
    title = "Key adenosine-related genes in OVARY DepMap (CERES)",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    # show x label once is already true (only 3 genes); keep it clean
    axis.text.x = element_text(face = "bold"),
    # avoid repeating long y labels across facets: keep only the first facet’s y text
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "4_4_heatmap_depmap_OVARY_ADA_ADORA2A_ADORA3.pdf"),
       p_markers, width = 9.0, height = 8.5)
