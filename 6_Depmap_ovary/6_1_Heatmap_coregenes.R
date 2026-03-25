#!/usr/bin/env Rscript
#load("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/5_3_HM_DepMap.RData")

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

DGE_PATH <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv"
TF_PATH  <- "~/Ovary_signatures/5_MRA/core_mra.tsv"

# ============================================================
# 1) Load consensus genes
# ============================================================

DGE_list_shared <- vroom(DGE_PATH, show_col_types = FALSE)
colnames(DGE_list_shared)
#[1] "gene"                "ensembl"             "AE"                 
#[4] "GTEx"                "GEO"                 "gtex"               
#[7] "ae"                  "geo"                 "n_sources"          
#[10] "sources"             "n_up"                "n_down"             
#[13] "consensus_direction" "log2FC_raw"          "log2FC_ae"          
#[16] "mean_log2FC"


cons_genes <- DGE_list_shared %>%
  filter(n_sources == 3,
         !is.na(consensus_direction),
         consensus_direction != "none") %>%
  transmute(gene_u = toupper(gene)) %>%
  distinct()

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
    dep_scaled_m11 = dep_raw ) %>% 
  #      2 * ((dep_raw - min(dep_raw, na.rm = TRUE)) /
  #             (max(dep_raw, na.rm = TRUE) - min(dep_raw, na.rm = TRUE))) - 1
  #  ) %>%
  ungroup()
colnames(crispr_ovary_all)
#[1] "depmap_id"      "gene"           "dependency"     "entrez_id"     
#[5] "gene_name"      "cell_line"      "gene_u"         "dep_raw"       
#[9] "dep_scaled_m11"

cat("Rows OVARY:", nrow(crispr_ovary_all), "\n")
#Rows OVARY: 1008388 

cat("Cell lines:", n_distinct(crispr_ovary_all$cell_line), "\n")
#Cell lines: 58 

# ============================================================
# 3) Build gene dependency matrix
# ============================================================

mat <- crispr_ovary_all %>%
  inner_join(cons_genes, by="gene_u") %>%
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
#Matrix dims: 611 x 58 

# ============================================================
# Heatmap + labels con flechas (sin traslape)
# ============================================================

# Ordenar matriz
mat_top <- mat[order(rowMeans(mat)), ]

# Escala de color
min_dep <- min(crispr_ovary_all$dep_raw, na.rm = TRUE)
max_dep <- max(crispr_ovary_all$dep_raw, na.rm = TRUE)

col_fun <- colorRamp2(
  c(min_dep, 0, max_dep),
  c("#B40426", "white", "#08306B")
)

# -------------------------------
# Dendrograma
# -------------------------------
row_dend <- as.dendrogram(hclust(dist(mat_top)))
row_dend <- reorder(row_dend, wts = -rowMeans(mat_top))
row_dend <- rev(row_dend)

# -------------------------------
# Funciones
# -------------------------------
get_leaves <- function(d) {
  if (is.leaf(d)) attr(d, "label")
  else unlist(lapply(d, get_leaves))
}

find_most_negative_clade <- function(dend, mat) {
  best_mean <- Inf
  best_genes <- NULL
  
  traverse <- function(node) {
    genes <- get_leaves(node)
    if (length(genes) > 1) {
      m <- mean(rowMeans(mat[genes, , drop=FALSE]), na.rm=TRUE)
      if (m < best_mean) {
        best_mean <<- m
        best_genes <<- genes
      }
    }
    if (!is.leaf(node)) lapply(node, traverse)
  }
  
  traverse(dend)
  best_genes
}

# Clado más negativo
genes_core <- find_most_negative_clade(row_dend, mat_top)

# Subir un nivel
get_parent_clade <- function(dend, target) {
  find_parent <- function(node) {
    if (!is.leaf(node)) {
      leaves <- get_leaves(node)
      if (all(target %in% leaves)) return(leaves)
      res <- lapply(node, find_parent)
      res <- res[!sapply(res, is.null)]
      if (length(res) > 0) return(res[[1]])
    }
    NULL
  }
  find_parent(dend)
}

genes_cluster <- get_parent_clade(row_dend, genes_core)

# -------------------------------
# TOP 10 genes más esenciales
# -------------------------------
gene_means <- rowMeans(mat_top)

top10_genes <- genes_cluster[
  order(gene_means[genes_cluster])[1:19]
]

# -------------------------------
# POSICIONES PARA anno_mark
# -------------------------------
row_index <- which(rownames(mat_top) %in% top10_genes)

# -------------------------------
# Heatmap + labels tipo flechas
# -------------------------------
ht <- Heatmap(
  mat_top,
  name = "Dependency_score",
  col = col_fun,
  cluster_rows = row_dend,
  cluster_columns = TRUE,
  
  show_row_names = FALSE,
  
  right_annotation = rowAnnotation(
    mark = anno_mark(
      at = row_index,
      labels = rownames(mat_top)[row_index],
      labels_gp = gpar(fontsize = 8, fontface = "bold"),
      link_width = unit(5, "mm"),
      link_gp = gpar(lwd = 0.8)
    )
  ),
  
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8,fontface = "bold"),
  column_title = "Gene dependency across ovarian cancer cell lines",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  
  rect_gp = gpar(col = NA),
  border = FALSE
)

# Dibujar
grid.newpage()
draw(ht, heatmap_legend_side = "right")

# Guardar
pdf("~/Ovary_signatures/6_Depmap_ovary/6_1_1_Heatmap_coregenes.pdf", 8.5, 8.5)
grid.newpage(); draw(ht, heatmap_legend_side="right"); dev.off()

png("~/Ovary_signatures/6_Depmap_ovary/6_1_1_Heatmap_coregenes.png",
    width=8.5, height=8.5, units="in", res=600)
grid.newpage(); draw(ht, heatmap_legend_side="right"); dev.off()



#save.image("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/5_3_HM_DepMap.RData")
