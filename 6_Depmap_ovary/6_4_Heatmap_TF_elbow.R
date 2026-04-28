#!/usr/bin/env Rscript
# ============================================================
# CRISPR dependency analysis – ovarian cancer
# ============================================================
#BiocManager::install("UCLouvain-CBIO/depmap")

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(effsize)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
})
options(stringsAsFactors = FALSE)

# ============================================================
# 1) Load regulones TF
# ============================================================

# All names of TF in regulones
Top_MR_DepMap <- readRDS("/STORAGE/csbig/jruiz/Ovary_data/6_MRA/6_0_1_Top_elbow_regulones_nes.rds")
MR_DGE <- Top_MR_DepMap %>% filter(isCore ==TRUE) %>% pull(TF) #11

Top_MR_DepMap <- Top_MR_DepMap %>% pull(TF)
length(Top_MR_DepMap) #27

# ============================================================
# 2) DepMap CRISPR ovarian cell lines
# ============================================================

crispr_ovary_all <- depmap_crispr() %>%
  filter(grepl("_OVARY$", cell_line)) %>%
  mutate(
    gene_u = toupper(gene_name),
    dep_raw = dependency
  ) %>%
  group_by(cell_line) %>%
  mutate(
    dep_scaled_m11 = dep_raw
  ) %>%
  ungroup()

colnames(crispr_ovary_all)

cat("Rows OVARY:", nrow(crispr_ovary_all), "\n")
cat("Cell lines:", n_distinct(crispr_ovary_all$cell_line), "\n")


# ============================================================
# 3) Build gene dependency matrix
# ============================================================

mat <- crispr_ovary_all %>%
  filter(gene_u %in% Top_MR_DepMap) %>%
  transmute(
    gene = gene_u,
    cell = sub("_OVARY$", "", cell_line),
    dep_scaled_m11
  ) %>%
  distinct() %>%
  pivot_wider(names_from = cell, values_from = dep_scaled_m11) %>%
  column_to_rownames("gene") %>%
  as.matrix()

mat[is.na(mat)] <- NA

cat("Matrix dims:", nrow(mat), "x", ncol(mat), "\n")


# ============================================================
# Heatmap 
# ============================================================

# -------------------------------
# ORDENAR MATRIZ (por dependencia media)
# -------------------------------
gene_means <- rowMeans(mat, na.rm = TRUE)
mat_top <- mat[order(gene_means, decreasing = FALSE), ]

# -------------------------------
# Scale
# -------------------------------
min_dep <- min(crispr_ovary_all$dep_raw, na.rm = TRUE)
max_dep <- max(crispr_ovary_all$dep_raw, na.rm = TRUE)

col_fun <- colorRamp2(
  c(min_dep, 0, max_dep),
  c("#B40426", "white", "#08306B")
)

# -------------------------------
# 🔥 FORMATO CONDICIONAL DE GENES
# -------------------------------
genes_ordered <- rownames(mat_top)

is_mr_dge <- genes_ordered %in% MR_DGE

row_names_gp <- gpar(
  fontsize = 9,
  fontface = ifelse(is_mr_dge, "bold", "plain")
)

# -------------------------------
# HEATMAP (TODOS LOS GENES)
# -------------------------------
ht <- Heatmap(
  mat_top,
  name = "Dependency_score",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  row_names_gp = row_names_gp,   
  
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  
  column_title = "Gene dependency across ovarian cancer cell lines",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  
  rect_gp = gpar(col = NA),
  border = FALSE
)

grid.newpage()
draw(ht, heatmap_legend_side = "right")


# ============================================================
# Save
# ============================================================

pdf("~/Ovary_signatures/6_Depmap_ovary/6_4_1_Heatmap_TF.pdf", 8.5, 8.5)
grid.newpage(); draw(ht, heatmap_legend_side="right"); dev.off()

png("~/Ovary_signatures/6_Depmap_ovary/6_4_1_Heatmap_TF.png",
    width=8.5, height=8.5, units="in", res=600)
grid.newpage(); draw(ht, heatmap_legend_side="right"); dev.off()

