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
  library(tibble)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Goal
# - Filter DepMap CRISPR dependency table to OVARY cell lines
# - Normalize "importance" PER cell line to [-1, 1] using ALL genes (pre-filter)
#     * +1 = most important (most essential) gene in that cell line
#     * -1 = least important (least essential) gene in that cell line
# - Keep only genes in DGE_list_shared with consensus_direction != "none"
# - Add consensus_direction and ovary_group columns
# - Save ONE output TSV
# - Plot heatmap using normalized importance [-1, 1]
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
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  "SLC28A1", "SLC29A1",
  "NT5E", "ENTPD1", "DPP4", "ADK", "ADA", "ADA2"
)))

# -------------------- OVARY line stratification (rule-of-thumb) -------------------- #
hgsoc_like <- c(
  "KURAMOCHI_OVARY", "OVCAR3_OVARY", "OVCAR4_OVARY", "OVCAR5_OVARY", "OVCAR8_OVARY",
  "OV90_OVARY", "CAOV3_OVARY", "CAOV4_OVARY", "OVCA420_OVARY", "COV362_OVARY",
  "NIHOVCAR3_OVARY"
)

non_epi_rare_germ <- c(
  "PA1_OVARY",
  "COV434_OVARY",
  "SCCOHT1_OVARY",
  "BIN67_OVARY"
)

assign_ovary_group <- function(cell_line) {
  if (cell_line %in% hgsoc_like) return("HGSOC-like")
  if (cell_line %in% non_epi_rare_germ) return("Non-epithelial/rare/germline")
  return("Other epithelial/mixed")
}

# -------------------- Load consensus table -------------------- #
DGE_list_shared <- vroom(DGE_PATH, show_col_types = FALSE)

cons_genes <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(!is.na(consensus_direction), consensus_direction != "none") %>%
  distinct(gene_u, consensus_direction)

cat("Consensus genes kept (direction != none): ", nrow(cons_genes), "\n", sep = "")

# -------------------- DepMap CRISPR (OVARY only) -------------------- #
# 1) Load OVARY subset (ALL genes)
# 2) Convert dependency -> importance (higher = more essential) by flipping sign
# 3) Normalize importance to [-1, 1] PER cell_line using ALL genes (pre-filter)
# 4) THEN filter to consensus genes

crispr_ovary_all <- depmap_crispr() %>%
  filter(grepl("_OVARY$", cell_line)) %>%
  mutate(
    gene_u = toupper(gene_name),
    ovary_group = vapply(cell_line, assign_ovary_group, character(1)),
    # DepMap dependency: more negative => more essential
    # Convert to "importance": higher => more essential
    importance_raw = -dependency
  )

cat("Rows in OVARY CRISPR (ALL genes): ", nrow(crispr_ovary_all), "\n", sep = "")
cat("Unique OVARY cell lines (ALL genes): ", dplyr::n_distinct(crispr_ovary_all$cell_line), "\n", sep = "")

# ---- Min-max normalize importance_raw to [0,1] within each cell line, then map to [-1,1]
# +1 = max importance (most essential), -1 = min importance (least essential)
crispr_ovary_all <- crispr_ovary_all %>%
  group_by(cell_line) %>%
  mutate(
    imp_min = min(importance_raw, na.rm = TRUE),
    imp_max = max(importance_raw, na.rm = TRUE),
    importance_01 = dplyr::if_else(
      is.finite(imp_min) & is.finite(imp_max) & (imp_max > imp_min),
      (importance_raw - imp_min) / (imp_max - imp_min),
      0.5  # if constant/degenerate, put in middle
    ),
    importance_m11 = 2 * importance_01 - 1
  ) %>%
  ungroup() %>%
  select(-imp_min, -imp_max)

# -------------------- Now filter to consensus genes (AFTER normalization) -------------------- #
crispr_ovary <- crispr_ovary_all %>%
  inner_join(cons_genes, by = "gene_u") %>%
  transmute(
    depmap_id,
    cell_line,
    ovary_group,
    gene = gene_u,
    dependency,            # original DepMap score
    importance_m11,        # NEW: normalized [-1,1], +1 most important per cell line
    consensus_direction
  )

cat("Unique genes in filtered DepMap table: ", length(unique(crispr_ovary$gene)), "\n", sep = "")
cat("Unique OVARY cell lines: ", length(unique(crispr_ovary$cell_line)), "\n", sep = "")

# Quick sanity counts per group
cellline_groups <- crispr_ovary %>%
  dplyr::select(cell_line, ovary_group) %>%
  dplyr::distinct() %>%
  tibble::as_tibble()

cat("\nCell lines per ovary_group:\n")
print(dplyr::count(cellline_groups, ovary_group, sort = TRUE))

# -------------------- Save output TSV -------------------- #
write_tsv(crispr_ovary, out_crispr)
cat("Saved: ", out_crispr, "\n", sep = "")

# ============================================================
# ===================== HEATMAP (normalized) =================
# ============================================================

# Remove "_OVARY" from cell line labels (for plotting)
crispr_plot <- crispr_ovary %>%
  mutate(cell = sub("_OVARY$", "", cell_line)) %>%
  select(gene, cell, ovary_group, importance_m11, consensus_direction)

# Build matrix genes x cells (use normalized importance)
mat_df <- crispr_plot %>%
  select(gene, cell, importance_m11) %>%
  distinct(gene, cell, .keep_all = TRUE) %>%
  pivot_wider(names_from = cell, values_from = importance_m11)

mat <- mat_df %>%
  as.data.frame()
rownames(mat) <- mat$gene
mat$gene <- NULL
mat <- as.matrix(mat)

# -------------------- Row annotation: consensus_direction -------------------- #
row_dir <- crispr_plot %>%
  distinct(gene, consensus_direction) %>%
  as.data.frame()
rownames(row_dir) <- row_dir$gene
row_dir$gene <- NULL
row_dir <- row_dir[rownames(mat), , drop = FALSE]

# -------------------- Column annotation: ovary_group -------------------- #
col_grp <- crispr_plot %>%
  distinct(cell, ovary_group) %>%
  as.data.frame()
rownames(col_grp) <- col_grp$cell
col_grp$cell <- NULL
col_grp <- col_grp[colnames(mat), , drop = FALSE]

# -------------------- Colors -------------------- #
# Since we normalized to [-1,1], use a fixed symmetric scale for comparability
col_fun <- circlize::colorRamp2(
  c(-1, 0, 1),
  c("#2C7BB6", "white", "#D7191C")
)

dir_cols <- c(up = "#1B9E77", down = "#D95F02")

# column group colors (simple, readable)
grp_levels <- unique(col_grp$ovary_group)
grp_cols <- setNames(
  c("#4E79A7", "#F28E2B", "#59A14F")[seq_along(grp_levels)],
  grp_levels
)

# Only label highlighted genes that exist
genes_present <- rownames(mat)
highlight_present <- intersect(highlight_genes, genes_present)
row_labels <- ifelse(genes_present %in% highlight_present, genes_present, "")

# -------------------- Annotations -------------------- #
ra <- rowAnnotation(
  direction = row_dir$consensus_direction,
  col = list(direction = dir_cols),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 9)
)

ca <- HeatmapAnnotation(
  ovary_group = col_grp$ovary_group,
  col = list(ovary_group = grp_cols),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9)
)

# -------------------- Heatmap -------------------- #
# Improve readability:
# - split columns by ovary_group
# - keep column clustering within groups (cluster_column_slices = TRUE)
# - hide row dendrogram (too many rows)
ht <- Heatmap(
  mat,
  name = "Importance\n(-1 to 1)",
  col = col_fun,
  na_col = "grey90",
  left_annotation = ra,
  top_annotation = ca,
  column_split = col_grp$ovary_group,
  cluster_column_slices = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  row_labels = row_labels,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(6, "cm"),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  column_title = "DepMap CRISPR normalized gene importance per OVARY cell line ([-1,1], +1 most essential)\nGenes from consensus DGE (direction != none)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    at = c(-1, -0.5, 0, 0.5, 1)
  )
)

# -------------------- Save PDF -------------------- #
pdf(out_pdf, width = 13, height = 14, useDingbats = FALSE)
grid.newpage()
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Saved: ", out_pdf, "\n", sep = "")

# -------------------- Save PNG -------------------- #
png(out_png, width = 3900, height = 4200, res = 300)
grid.newpage()
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("Saved: ", out_png, "\n", sep = "")

cat("Highlighted genes present (labeled): ", paste(highlight_present, collapse = ", "), "\n", sep = "")
cat("Done.\n")





#####################
# ============================================================
# ========== EXTRA HEATMAPS: Master Regulators / Targets ======
# Requirements:
# 1) Heatmap for TOP 10 master regulators (TF) by rank_by_interest
#    - keep only rows with non-NA core_signature_score
#    - plot those TF genes in DepMap (same style as main heatmap)
# 2) Heatmap for ALL unique interest_target_genes across the table
#    - plot those genes in DepMap (same style)
# Notes:
# - Both use importance_m11 (normalized [-1,1]) already computed pre-filter.
# - Both keep the same column ordering + annotations (ovary_group split).
# ============================================================

# -------------------- Load MR table -------------------- #
masters_reg_path <- "~/Ovary_signatures/2_DEG_GEO/core_mrs_interest_annot.tsv"
masters_reg <- vroom::vroom(masters_reg_path, show_col_types = FALSE)

# -------------------- Helper: build and save heatmap from a gene set -------------------- #
make_heatmap_for_genes <- function(genes_vec,
                                   title,
                                   out_pdf_path,
                                   out_png_path,
                                   crispr_ovary_all_obj,
                                   cellline_groups_df,
                                   highlight_genes_vec = character(0)) {
  genes_vec <- unique(toupper(genes_vec))
  genes_vec <- genes_vec[!is.na(genes_vec) & genes_vec != ""]
  
  # subset from ALL (already normalized)
  df <- crispr_ovary_all_obj %>%
    mutate(gene_u = toupper(gene_u)) %>%
    filter(gene_u %in% genes_vec) %>%
    transmute(
      cell_line,
      gene = gene_u,
      importance_m11
    )
  
  if (nrow(df) == 0) {
    warning("No matching genes found in DepMap for: ", title)
    return(invisible(NULL))
  }
  
  # Keep only OVARY cell lines; remove suffix for plotting
  df_plot <- df %>%
    mutate(cell = sub("_OVARY$", "", cell_line)) %>%
    select(gene, cell, importance_m11)
  
  # matrix
  mat_df2 <- df_plot %>%
    distinct(gene, cell, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = cell, values_from = importance_m11)
  
  mat2 <- mat_df2 %>% as.data.frame()
  rownames(mat2) <- mat2$gene
  mat2$gene <- NULL
  mat2 <- as.matrix(mat2)
  
  # Ensure columns order matches your existing main heatmap (same as mat)
  # If 'mat' exists in env (from main heatmap), use it; otherwise use current order.
  if (exists("mat") && is.matrix(get("mat"))) {
    main_cols <- colnames(get("mat"))
    common_cols <- intersect(main_cols, colnames(mat2))
    mat2 <- mat2[, common_cols, drop = FALSE]
  }
  
  # Column annotation from cellline_groups_df (cell_line has _OVARY, mat has no suffix)
  col_grp2 <- cellline_groups_df %>%
    mutate(cell = sub("_OVARY$", "", cell_line)) %>%
    select(cell, ovary_group) %>%
    distinct() %>%
    as.data.frame()
  rownames(col_grp2) <- col_grp2$cell
  col_grp2$cell <- NULL
  col_grp2 <- col_grp2[colnames(mat2), , drop = FALSE]
  
  # colors: fixed [-1,1]
  col_fun2 <- circlize::colorRamp2(c(-1, 0, 1), c("#2C7BB6", "white", "#D7191C"))
  
  # group colors (re-use if exist, else define)
  if (!exists("grp_cols")) {
    grp_levels2 <- unique(col_grp2$ovary_group)
    grp_cols2 <- setNames(
      c("#4E79A7", "#F28E2B", "#59A14F")[seq_along(grp_levels2)],
      grp_levels2
    )
  } else {
    grp_cols2 <- get("grp_cols")
  }
  
  ca2 <- ComplexHeatmap::HeatmapAnnotation(
    ovary_group = col_grp2$ovary_group,
    col = list(ovary_group = grp_cols2),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9)
  )
  
  # row labels: show all (small set) OR only highlight if you want
  genes_present2 <- rownames(mat2)
  if (length(highlight_genes_vec) > 0) {
    hl2 <- intersect(toupper(highlight_genes_vec), genes_present2)
    row_labels2 <- ifelse(genes_present2 %in% hl2, genes_present2, "")
    show_row_names2 <- TRUE
    row_names_gp2 <- gpar(fontsize = 9)
  } else {
    row_labels2 <- genes_present2
    show_row_names2 <- TRUE
    row_names_gp2 <- gpar(fontsize = 9)
  }
  
  ht2 <- ComplexHeatmap::Heatmap(
    mat2,
    name = "Importance\n(-1 to 1)",
    col = col_fun2,
    na_col = "grey90",
    top_annotation = ca2,
    column_split = col_grp2$ovary_group,
    cluster_column_slices = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_labels = row_labels2,
    show_row_names = show_row_names2,
    row_names_gp = row_names_gp2,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8),
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      at = c(-1, -0.5, 0, 0.5, 1)
    )
  )
  
  # save
  pdf(out_pdf_path, width = 12, height = max(6, 0.25 * nrow(mat2) + 4), useDingbats = FALSE)
  grid.newpage()
  draw(ht2, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("Saved: ", out_pdf_path, "\n", sep = "")
  
  png(out_png_path, width = 3600, height = 1800 + 120 * nrow(mat2), res = 300)
  grid.newpage()
  draw(ht2, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("Saved: ", out_png_path, "\n", sep = "")
  
  invisible(list(mat = mat2, ht = ht2))
}

# -------------------- 1) TOP 10 master regulators (TF) -------------------- #
top10_mr <- masters_reg %>%
  filter(!is.na(core_signature_score)) %>%              # your rule
  arrange(rank_by_interest) %>%
  distinct(TF, .keep_all = TRUE) %>%
  slice_head(n = 10)

top10_tfs <- unique(toupper(top10_mr$TF))

cat("\nTop 10 MRs (TF) used for heatmap:\n")
print(top10_mr %>% select(TF, rank_by_interest, core_signature_score, TF_expression_status))

out_pdf_mr10 <- file.path(OUT_DIR, "5_3_DepMap_OVARY_heatmap_top10_master_regulators.pdf")
out_png_mr10 <- file.path(OUT_DIR, "5_3_DepMap_OVARY_heatmap_top10_master_regulators.png")

make_heatmap_for_genes(
  genes_vec = top10_tfs,
  title = "DepMap OVARY — Top 10 Master Regulators (TF) by rank_by_interest (core_signature_score not NA)\nImportance normalized [-1,1] (+1 most essential per cell line)",
  out_pdf_path = out_pdf_mr10,
  out_png_path = out_png_mr10,
  crispr_ovary_all_obj = crispr_ovary_all,
  cellline_groups_df = cellline_groups
)

# -------------------- 2) ALL unique interest_target_genes -------------------- #
# split on ';' and clean
all_targets <- masters_reg %>%
  filter(!is.na(interes_target_genes), interes_target_genes != "") %>%
  pull(interes_target_genes) %>%
  strsplit(";", fixed = TRUE) %>%
  unlist(use.names = FALSE)

all_targets <- unique(toupper(trimws(all_targets)))
all_targets <- all_targets[!is.na(all_targets) & all_targets != ""]

cat("\nUnique interest_target_genes (all rows): ", length(all_targets), "\n", sep = "")

out_pdf_targets <- file.path(OUT_DIR, "5_4_DepMap_OVARY_heatmap_all_interest_target_genes.pdf")
out_png_targets <- file.path(OUT_DIR, "5_4_DepMap_OVARY_heatmap_all_interest_target_genes.png")

make_heatmap_for_genes(
  genes_vec = all_targets,
  title = paste0(
    "DepMap OVARY — All unique interest_target_genes (n=", length(all_targets), ")\n",
    "Importance normalized [-1,1] (+1 most essential per cell line)"
  ),
  out_pdf_path = out_pdf_targets,
  out_png_path = out_png_targets,
  crispr_ovary_all_obj = crispr_ovary_all,
  cellline_groups_df = cellline_groups,
  highlight_genes_vec = highlight_genes  # will label these if present (optional)
)

cat("Done: extra MR/target heatmaps.\n")




