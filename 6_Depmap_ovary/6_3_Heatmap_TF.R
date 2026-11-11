#!/usr/bin/env Rscript
BASE <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))  # render Unicode minus signs correctly

# ============================================================
# CRISPR dependency analysis – ovarian cancer
# ============================================================
#BiocManager::install("UCLouvain-CBIO/depmap")

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
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
# 1) Load consensus genes
# ============================================================

# Names of all TFs in regulons
Top_MR_DepMap <- readRDS(paste0(BASE, "/6_Depmap_ovary/6_2_1_regulons_sig_filtered.rds"))
Top_MR_DepMap <- names(Top_MR_DepMap)
length(Top_MR_DepMap) #46


MR_DGE <- readRDS(file = paste0(BASE, "/6_Depmap_ovary/6_2_0_names_TF_DEG_TopNES.rds"))
length(MR_DGE) #23
# 
MR_DGE <- toupper(MR_DGE)

# ============================================================
# 2) DepMap CRISPR ovarian cell lines
# ============================================================
# DepMap release 24Q4 (CRISPRGeneEffect.csv = Chronos gene effect; Model.csv), downloaded from
# depmap.org. The Bioconductor 'depmap' package used previously ships release 22Q2.
DEPMAP_DIR <- Sys.getenv("DEPMAP_DIR", paste0(BASE, "/data/depmap_24Q4"))
models <- vroom(file.path(DEPMAP_DIR, "Model.csv"), show_col_types = FALSE)
ovary_lines <- models %>%
  filter(OncotreePrimaryDisease == "Ovarian Epithelial Tumor") %>%
  dplyr::select(depmap_id = ModelID, cell_line = StrippedCellLineName)

gene_effect <- vroom(file.path(DEPMAP_DIR, "CRISPRGeneEffect.csv"), show_col_types = FALSE)
names(gene_effect)[1] <- "depmap_id"
crispr_ovary_all <- gene_effect %>%
  inner_join(ovary_lines, by = "depmap_id") %>%
  pivot_longer(-c(depmap_id, cell_line), names_to = "gene_name", values_to = "dependency") %>%
  mutate(
    gene_name = sub(" \\(\\d+\\)$", "", gene_name),
    gene_u = toupper(gene_name),
    dep_raw = dependency,
    dep_scaled_m11 = dep_raw
  )

cat("Rows OVARY:", nrow(crispr_ovary_all), "\n")
cat("Cell lines:", n_distinct(crispr_ovary_all$cell_line), "\n")

colnames(crispr_ovary_all)


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

# Drop regulators without a 24Q4 gene effect in most lines (HMGB3: multi-aligning guides removed)
na_frac <- rowMeans(is.na(mat))
print(sort(na_frac[na_frac > 0], decreasing = TRUE))
mat <- mat[na_frac <= 0.5, , drop = FALSE]

cat("Matrix dims:", nrow(mat), "x", ncol(mat), "\n")


# ============================================================
# Heatmap 
# ============================================================

# -------------------------------
# ORDER MATRIX (by mean dependency)
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
# Conditional row formatting
# -------------------------------
genes_ordered <- rownames(mat_top)

is_mr_dge <- genes_ordered %in% MR_DGE

row_names_gp <- gpar(
  fontsize = 9,
  fontface = ifelse(is_mr_dge, "bold", "plain")
)

# -------------------------------
# HEATMAP (ALL GENES)
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
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  
  rect_gp = gpar(col = NA),
  border = FALSE
)

grid.newpage()
draw(ht, heatmap_legend_side = "right")


# ============================================================
# Save
# ============================================================

pdf(paste0(BASE, "/6_Depmap_ovary/6_3_1_Heatmap_TF.pdf"), 8.5, 8.5)
grid.newpage(); draw(ht, heatmap_legend_side="right"); dev.off()

png(paste0(BASE, "/6_Depmap_ovary/6_3_1_Heatmap_TF.png"),
    width=8.5, height=8.5, units="in", res=600)
grid.newpage(); draw(ht, heatmap_legend_side="right"); dev.off()




resultado <- data.frame(gene = rownames(mat), media = rowMeans(mat, na.rm = TRUE), 
                        mediana = apply(mat, 1, median, na.rm = TRUE)) %>% arrange(desc(media))

median <- resultado %>% filter(mediana <= -0.5) %>% pull(gene)






