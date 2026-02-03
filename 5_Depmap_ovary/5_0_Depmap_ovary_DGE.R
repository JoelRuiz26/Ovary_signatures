#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(readr)      # <- necesario para write_tsv()
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Goal
# - DepMap CRISPR OVARY: normalize importance PER cell line to [-1,1] using ALL genes (pre-filter)
# - Main heatmap: ONLY first 50 DGE genes + highlight_genes
# - TF heatmap: ALL TFs from your image (hard-coded, no ranking)
# - Interest heatmap: ONLY highlight_genes (do NOT extract from masters table)
# ============================================================

# -------------------- Paths -------------------- #
DGE_PATH <- "/home/jruiz/Ovary_signatures/3_Consensus_DGE_analysis/3_0_1_genes_concordant.tsv"
OUT_DIR  <- "/home/jruiz/Ovary_signatures/5_Depmap_ovary"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out_crispr <- file.path(OUT_DIR, "5_1_depmap_CRISPR_OVARY_genes_in_DGE_list_shared.tsv")

out_pdf_main <- file.path(OUT_DIR, "5_2_DepMap_OVARY_heatmap_DGE_top50_plus_highlight.pdf")
out_png_main <- file.path(OUT_DIR, "5_2_DepMap_OVARY_heatmap_DGE_top50_plus_highlight.png")

out_pdf_tf   <- file.path(OUT_DIR, "5_3_DepMap_OVARY_heatmap_TF_ALL_from_image.pdf")
out_png_tf   <- file.path(OUT_DIR, "5_3_DepMap_OVARY_heatmap_TF_ALL_from_image.png")

out_pdf_int  <- file.path(OUT_DIR, "5_4_DepMap_OVARY_heatmap_interest_highlight_genes.pdf")
out_png_int  <- file.path(OUT_DIR, "5_4_DepMap_OVARY_heatmap_interest_highlight_genes.png")

# -------------------- Genes to highlight -------------------- #
highlight_genes <- unique(toupper(c(
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  "SLC28A1", "SLC29A1",
  "NT5E", "ENTPD1", "DPP4", "ADK", "ADA", "ADA2"
)))

# -------------------- OVARY line stratification -------------------- #
hgsoc_like <- c(
  "KURAMOCHI_OVARY", "OVCAR3_OVARY", "OVCAR4_OVARY", "OVCAR5_OVARY", "OVCAR8_OVARY",
  "OV90_OVARY", "CAOV3_OVARY", "CAOV4_OVARY", "OVCA420_OVARY", "COV362_OVARY",
  "NIHOVCAR3_OVARY"
)

non_epi_rare_germ <- c("PA1_OVARY", "COV434_OVARY", "SCCOHT1_OVARY", "BIN67_OVARY")

assign_ovary_group <- function(cell_line) {
  if (cell_line %in% hgsoc_like) return("HGSOC-like")
  if (cell_line %in% non_epi_rare_germ) return("Non-epithelial/rare/germline")
  return("Other epithelial/mixed")
}

# ============================================================
# Load consensus DGE table
# ============================================================
DGE_list_shared <- vroom(DGE_PATH, show_col_types = FALSE)

cons_genes <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(!is.na(consensus_direction), consensus_direction != "none") %>%
  distinct(gene_u, consensus_direction)

cat("Consensus genes kept (direction != none): ", nrow(cons_genes), "\n", sep = "")

# ============================================================
# Select TOP 25 UP + TOP 25 DOWN DGE by effect size (|logFC|)
# ============================================================
top25_up <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(consensus_direction == "up", !is.na(logFC)) %>%
  arrange(desc(abs(logFC))) %>%
  distinct(gene_u, .keep_all = TRUE) %>%
  slice_head(n = 25) %>%
  pull(gene_u)

top25_down <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(consensus_direction == "down", !is.na(logFC)) %>%
  arrange(desc(abs(logFC))) %>%
  distinct(gene_u, .keep_all = TRUE) %>%
  slice_head(n = 25) %>%
  pull(gene_u)

top50_dge <- unique(c(top25_up, top25_down))
genes_main <- unique(c(top50_dge, highlight_genes))

cat("Top DGE selected:\n",
    " - UP:", length(top25_up), "\n",
    " - DOWN:", length(top25_down), "\n",
    " - TOTAL:", length(top50_dge), "\n", sep = "")
cat("Main heatmap genes (top50 + highlight): ", length(genes_main), "\n", sep = "")

# ============================================================
# DepMap CRISPR (OVARY) + normalization to [-1,1] per cell line
# ============================================================
crispr_ovary_all <- depmap_crispr() %>%
  filter(grepl("_OVARY$", cell_line)) %>%
  mutate(
    gene_u = toupper(gene_name),
    ovary_group = vapply(cell_line, assign_ovary_group, character(1)),
    importance_raw = -dependency
  )

cat("Rows in OVARY CRISPR (ALL genes): ", nrow(crispr_ovary_all), "\n", sep = "")
cat("Unique OVARY cell lines (ALL genes): ", n_distinct(crispr_ovary_all$cell_line), "\n", sep = "")

crispr_ovary_all <- crispr_ovary_all %>%
  group_by(cell_line) %>%
  mutate(
    imp_min = min(importance_raw, na.rm = TRUE),
    imp_max = max(importance_raw, na.rm = TRUE),
    importance_01 = if_else(
      is.finite(imp_min) & is.finite(imp_max) & (imp_max > imp_min),
      (importance_raw - imp_min) / (imp_max - imp_min),
      0.5
    ),
    importance_m11 = 2 * importance_01 - 1
  ) %>%
  ungroup() %>%
  dplyr::select(-imp_min, -imp_max)

# ============================================================
# Save TSV (full filtered to consensus genes, as before)
# ============================================================
crispr_ovary <- crispr_ovary_all %>%
  inner_join(cons_genes, by = "gene_u") %>%
  transmute(
    depmap_id,
    cell_line,
    ovary_group,
    gene = gene_u,
    dependency,
    importance_m11,
    consensus_direction
  )

write_tsv(crispr_ovary, out_crispr)
cat("Saved: ", out_crispr, "\n", sep = "")

cellline_groups <- crispr_ovary_all %>%
  distinct(cell_line, ovary_group) %>%
  as_tibble()

cat("\nCell lines per ovary_group:\n")
print(cellline_groups %>% count(ovary_group, sort = TRUE))

# ============================================================
# Plot helpers
# ============================================================
open_png_device <- function(path, width_px, height_px, res = 300) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width_px, height = height_px, units = "px", res = res, scaling = 1)
    return(invisible("ragg"))
  }
  png(filename = path, width = width_px, height = height_px, res = res, type = "cairo")
  invisible("cairo")
}

build_mat <- function(df_all, genes_vec) {
  genes_vec <- unique(toupper(genes_vec))
  df <- df_all %>%
    filter(gene_u %in% genes_vec) %>%
    transmute(cell = sub("_OVARY$", "", cell_line),
              gene = gene_u,
              importance_m11)
  
  mat_df <- df %>%
    distinct(gene, cell, .keep_all = TRUE) %>%
    pivot_wider(names_from = cell, values_from = importance_m11)
  
  mat <- as.data.frame(mat_df)
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  mat <- as.matrix(mat)
  
  if (anyNA(mat)) mat[is.na(mat)] <- 0
  mat
}

get_col_grp_for_mat <- function(mat, cellline_groups_df) {
  col_grp <- cellline_groups_df %>%
    mutate(cell = sub("_OVARY$", "", cell_line)) %>%
    dplyr::select(cell, ovary_group) %>%
    distinct() %>%
    as.data.frame()
  rownames(col_grp) <- col_grp$cell
  col_grp$cell <- NULL
  col_grp[colnames(mat), , drop = FALSE]
}

make_heatmap <- function(mat, col_grp, title, out_pdf, out_png, show_row_names = TRUE) {
  
  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2C7BB6", "white", "#D7191C"))
  
  grp_levels <- unique(col_grp$ovary_group)
  grp_cols <- setNames(c("#4E79A7", "#F28E2B", "#59A14F")[seq_along(grp_levels)], grp_levels)
  
  ca <- HeatmapAnnotation(
    ovary_group = col_grp$ovary_group,
    col = list(ovary_group = grp_cols),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
  )
  
  ht <- Heatmap(
    mat,
    name = "Importance\n(-1 to 1)",
    col = col_fun,
    top_annotation = ca,
    column_split = col_grp$ovary_group,
    cluster_column_slices = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    show_row_names = show_row_names,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8),
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1)),
    rect_gp = gpar(col = NA),
    border = FALSE,
    use_raster = FALSE
  )
  
  pdf_h <- min(10, max(6, 0.18 * nrow(mat) + 3))
  pdf(out_pdf, width = 12, height = pdf_h, useDingbats = FALSE)
  grid.newpage()
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("Saved: ", out_pdf, "\n", sep = "")
  
  height_px <- min(2600, max(1600, 900 + 45 * nrow(mat)))
  dev_used <- open_png_device(out_png, width_px = 3600, height_px = height_px, res = 300)
  grid.newpage()
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("Saved: ", out_png, " (device=", dev_used, ")\n", sep = "")
  
  invisible(ht)
}

# ============================================================
# 1) MAIN heatmap
# ============================================================
mat_main <- build_mat(crispr_ovary_all, genes_main)
col_grp_main <- get_col_grp_for_mat(mat_main, cellline_groups)

make_heatmap(
  mat = mat_main,
  col_grp = col_grp_main,
  title = paste0(
    "DepMap OVARY — DGE top50 (25 UP + 25 DOWN by |logFC|) + highlight genes (n=", nrow(mat_main), ")\n",
    "Importance normalized [-1,1] (+1 most essential per cell line)"
  ),
  out_pdf = out_pdf_main,
  out_png = out_png_main,
  show_row_names = TRUE
)

# ============================================================
# 2) TF heatmap: ALL TFs from your image (hard-coded order)
# ============================================================
tf_all_from_image <- toupper(c(
  "FLI1","ZEB2","PRDM8","PRDM10","FGF2","LYL1","ZNF624","SRSF2","HNRNPC","HDGF",
  "RELA","LDB2","MAZ","HSPA8","SFPQ","JUP","MEF2C","YY1AP1","SF1","PML",
  "HNRNPU","HNRNPK","ETS1","SNAI2","MUC1","CSDE1","GLI1","ZNF366","ETV1","PRDM1",
  "ZNF365","ZNF319","ATF7-NPFF1","XRN2","SATB2","ESR2","PICALM","CAMK4","FBXW7","ZMAT4"
))

mat_tf <- build_mat(crispr_ovary_all, tf_all_from_image)
col_grp_tf <- get_col_grp_for_mat(mat_tf, cellline_groups)

make_heatmap(
  mat = mat_tf,
  col_grp = col_grp_tf,
  title = paste0(
    "DepMap OVARY — TFs from image (n=", nrow(mat_tf), "; no re-ranking)\n",
    "Importance normalized [-1,1] (+1 most essential per cell line)"
  ),
  out_pdf = out_pdf_tf,
  out_png = out_png_tf,
  show_row_names = TRUE
)

# ============================================================
# 3) Interest heatmap: ONLY highlight_genes
# ============================================================
mat_int <- build_mat(crispr_ovary_all, highlight_genes)
col_grp_int <- get_col_grp_for_mat(mat_int, cellline_groups)

make_heatmap(
  mat = mat_int,
  col_grp = col_grp_int,
  title = paste0(
    "DepMap OVARY — Interest genes (highlight_genes) (n=", nrow(mat_int), ")\n",
    "Importance normalized [-1,1] (+1 most essential per cell line)"
  ),
  out_pdf = out_pdf_int,
  out_png = out_png_int,
  show_row_names = TRUE
)

cat("Done.\n")
