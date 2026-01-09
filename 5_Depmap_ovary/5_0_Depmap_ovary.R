#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(readr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Goal
# - Filter DepMap CRISPR dependency table to OVARY cell lines
# - Keep only genes in your DGE_list_shared with consensus_direction != "none"
# - Add consensus_direction (up/down) as a column to the DepMap table
# - Save ONE output TSV
# - Plot heatmap: 3391 genes (rows) x 58 OVARY cell lines (cols)
#   * No row gene names (too many)
#   * Only highlight selected genes with labels IF they exist
#   * Save PDF + PNG (professional)
# ============================================================

# -------------------- Paths -------------------- #
DGE_PATH <- "/home/jruiz/Ovary_signatures/3_Consensus_DGE_analysis/3_0_1_genes_concordant.tsv"
OUT_DIR  <- "/home/jruiz/Ovary_signatures/5_Depmap_ovary"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out_crispr <- file.path(OUT_DIR, "5_1_depmap_CRISPR_OVARY_genes_in_DGE_list_shared.tsv")

out_pdf <- file.path(OUT_DIR, "5_2_DepMap_OVARY_dependency_heatmap.pdf")
out_png <- file.path(OUT_DIR, "5_2_DepMap_OVARY_dependency_heatmap.png")

# -------------------- Genes to highlight (labels only if present) -------------------- #
highlight_genes <- unique(toupper(c(
  # Receptores
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  # Transportadores
  "SLC28A1", "SLC29A1",
  # Enzimas
  "NT5E", "ENTPD1", "DPP4", "ADK", "ADA", "ADA2"
)))

# -------------------- Load consensus table -------------------- #
DGE_list_shared <- vroom(DGE_PATH, show_col_types = FALSE)

cons_genes <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(!is.na(consensus_direction), consensus_direction != "none") %>%
  distinct(gene_u, consensus_direction)

cat("Consensus genes kept (direction != none): ", nrow(cons_genes), "\n", sep = "")

# -------------------- DepMap CRISPR (OVARY only) -------------------- #
# NOTE: This downloads/loads DepMap data via depmap package.
crispr_ovary <- depmap_crispr() %>%
  filter(grepl("_OVARY$", cell_line)) %>%
  mutate(gene_u = toupper(gene_name)) %>%
  inner_join(cons_genes, by = "gene_u") %>%
  transmute(
    depmap_id,
    cell_line,
    gene = gene_u,                 # keep uppercase for stable matching/labels
    dependency,
    consensus_direction
  )

cat("Unique genes in filtered DepMap table: ", length(unique(crispr_ovary$gene)), "\n", sep = "")
cat("Unique OVARY cell lines: ", length(unique(crispr_ovary$cell_line)), "\n", sep = "")

# -------------------- Save output TSV -------------------- #
write_tsv(crispr_ovary, out_crispr)
cat("Saved: ", out_crispr, "\n", sep = "")

# ===================== HEATMAP ===================== #

# Remove "_OVARY" from cell line labels (for plotting)
crispr_plot <- crispr_ovary %>%
  mutate(cell = sub("_OVARY$", "", cell_line)) %>%
  select(gene, cell, dependency, consensus_direction)

# Build matrix genes x cells
mat_df <- crispr_plot %>%
  select(gene, cell, dependency) %>%
  distinct(gene, cell, .keep_all = TRUE) %>%
  pivot_wider(names_from = cell, values_from = dependency)

mat <- mat_df %>%
  as.data.frame()

rownames(mat) <- mat$gene
mat$gene <- NULL
mat <- as.matrix(mat)

# Row annotation: consensus_direction (up/down) per gene
row_dir <- crispr_plot %>%
  distinct(gene, consensus_direction) %>%
  as.data.frame()
rownames(row_dir) <- row_dir$gene
row_dir$gene <- NULL
row_dir <- row_dir[rownames(mat), , drop = FALSE]

# Dynamic color range (robust to outliers):
# use 1% and 99% quantiles, then make it symmetric around 0 for a balanced diverging map
q <- as.numeric(stats::quantile(mat, probs = c(0.01, 0.99), na.rm = TRUE))
max_abs <- max(abs(q))
col_fun <- circlize::colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#2C7BB6", "white", "#D7191C")
)

# Only label highlighted genes that exist
genes_present <- rownames(mat)
highlight_present <- intersect(highlight_genes, genes_present)
row_labels <- ifelse(genes_present %in% highlight_present, genes_present, "")

# A small row annotation for direction (optional but scientific/useful)
dir_cols <- c(up = "#1B9E77", down = "#D95F02")
ra <- rowAnnotation(
  direction = row_dir$consensus_direction,
  col = list(direction = dir_cols),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 9)
)

# Heatmap object
ht <- Heatmap(
  mat,
  name = "CRISPR\ndependency",
  col = col_fun,
  na_col = "grey90",
  left_annotation = ra,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,          # cleaner with 3391 rows
  show_column_dend = TRUE,
  row_labels = row_labels,        # blank for most, names only for highlighted
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(5.5, "cm"),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  column_title = "DepMap CRISPR dependency (OVARY cell lines) — genes from consensus DGE (direction != none)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# -------------------- Save PDF -------------------- #
pdf(out_pdf, width = 12, height = 14, useDingbats = FALSE)
grid.newpage()
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Saved: ", out_pdf, "\n", sep = "")

# -------------------- Save PNG -------------------- #
png(out_png, width = 3600, height = 4200, res = 300)
grid.newpage()
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Saved: ", out_png, "\n", sep = "")

cat("Highlighted genes present (labeled): ", paste(highlight_present, collapse = ", "), "\n", sep = "")
cat("Done.\n")
