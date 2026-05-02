#!/usr/bin/env Rscript
#load("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/5_3_HM_DepMap.RData")

# ============================================================
# CRISPR dependency analysis – ovarian cancer
# ============================================================

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(patchwork)
})

options(stringsAsFactors = FALSE)


# ============================================================
# Paths
# ============================================================

DGE_PATH <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv"
OUT_DIR  <- "~/Ovary_signatures/6_Depmap_ovary/"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ============================================================
# 1) Load consensus genes
# ============================================================

DGE_list_shared <- vroom(DGE_PATH, show_col_types = FALSE)

cons_genes <- DGE_list_shared %>%
  filter(n_sources == 3,
         !is.na(consensus_direction),
         consensus_direction != "none") %>%
  transmute(gene_u = toupper(gene)) %>%
  distinct()

# ============================================================
# 2) Load DepMap CRISPR ovarian cell lines
# ============================================================
meta <- depmap_metadata()

ovary_lines <- meta %>%
  filter(lineage == "ovary") %>%
  filter(primary_disease== "Ovarian Cancer") %>% 
  pull(depmap_id)

crispr_ovary_all <- depmap_crispr() %>%
  filter(depmap_id %in% ovary_lines) %>%
  mutate(
    gene_u = toupper(gene_name),
    dep_raw = dependency
  )%>%
  group_by(cell_line) %>%
  mutate(dep_scaled_m11 = dep_raw) %>%
  ungroup()

cat("Rows OVARY:", nrow(crispr_ovary_all), "\n")
cat("Cell lines:", n_distinct(crispr_ovary_all$cell_line), "\n")

# ============================================================
# 3) Build gene dependency matrix
# ============================================================

mat <- crispr_ovary_all %>%
  inner_join(cons_genes, by = "gene_u") %>%
  transmute(
    gene = gene_u,
    cell = sub("_OVARY$", "", cell_line),
    dep_scaled_m11
  ) %>%
  distinct() %>%
  pivot_wider(names_from = cell, values_from = dep_scaled_m11) %>%
  column_to_rownames("gene") %>%
  as.matrix()

cat("Matrix dims:", nrow(mat), "x", ncol(mat), "\n")

# ============================================================
# 4) Heatmap (ordered by mean dependency)
# ============================================================

gene_means <- rowMeans(mat, na.rm = TRUE)
mat_top <- mat[order(gene_means, decreasing = FALSE), ]

min_dep <- min(crispr_ovary_all$dep_raw, na.rm = TRUE)
max_dep <- max(crispr_ovary_all$dep_raw, na.rm = TRUE)

col_fun <- colorRamp2(
  c(min_dep, 0, max_dep),
  c("#B40426", "white", "#08306B")
)

top_genes <- names(sort(rowMeans(mat_top, na.rm = TRUE)))[1:20]
row_index <- which(rownames(mat_top) %in% top_genes)

ht <- Heatmap(
  mat_top,
  name = "Dependency_score",
  col = col_fun,
  cluster_rows = FALSE,
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
  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_title = "Gene dependency across ovarian cancer cell lines",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  rect_gp = gpar(col = NA),
  border = FALSE
)

# Save heatmap
pdf(paste0(OUT_DIR, "6_1_1_Heatmap_coregenes.pdf"), 8.5, 8.5)
grid.newpage(); draw(ht, heatmap_legend_side = "right"); dev.off()

png(paste0(OUT_DIR, "6_1_1_Heatmap_coregenes.png"),
    width = 8.5, height = 8.5, units = "in", res = 600)
grid.newpage(); draw(ht, heatmap_legend_side = "right"); dev.off()

# ============================================================
# 5) Ranked dependency plot
# ============================================================

df_all_summary <- as.data.frame(mat) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "dep") %>%
  group_by(gene) %>%
  summarise(
    mean_dep  = mean(dep, na.rm = TRUE),
    sd_dep    = sd(dep, na.rm = TRUE),
    median_dep= median(dep, na.rm = TRUE),
    pct_neg   = mean(dep < -0.5, na.rm = TRUE),
    cv        = sd_dep / abs(mean_dep),
    n_obs     = sum(!is.na(dep)),
    .groups   = "drop"
  ) %>%
  arrange(mean_dep) %>%
  mutate(rank = row_number())

df_crisp_deg <- df_all_summary %>% left_join(DGE_list_shared , by = "gene") %>% 
  dplyr::select(gene, median_dep,mean_log2FC,consensus_direction) %>% 
  mutate(direction = consensus_direction,
         median_log2FC = mean_log2FC) %>% dplyr::select(-consensus_direction, -mean_log2FC)
vroom_write(df_crisp_deg, paste0(OUT_DIR, "6_1_0_Dependency_DEG_core.tsv"))


TOP_N <- 10
zoom_data <- df_all_summary %>% filter(rank <= TOP_N)

# Color scales (data-driven)
pal_rank_color <- scale_color_gradientn(
  colours = c("#B40426", "#FFDD57", "#2166AC"),
  values = scales::rescale(c(
    min(df_all_summary$mean_dep, na.rm = TRUE),
    0,
    max(df_all_summary$mean_dep, na.rm = TRUE)
  )),
  name = "Mean dep. score"
)

pal_rank_fill <- scale_fill_gradientn(
  colours = c("#B40426", "#FFDD57", "#2166AC"),
  values = scales::rescale(c(
    min(df_all_summary$mean_dep, na.rm = TRUE),
    0,
    max(df_all_summary$mean_dep, na.rm = TRUE)
  )),
  name = "Mean dep. score"
)

# Main plot
p_main <- ggplot(df_all_summary, aes(rank, mean_dep, color = mean_dep)) +
  geom_ribbon(aes(ymin = mean_dep - sd_dep, ymax = mean_dep + sd_dep, fill = mean_dep),
              alpha = 0.18, color = NA) +
  pal_rank_fill +
  geom_line(linewidth = 0.35, alpha = 0.45, color = "grey50") +
  geom_point(size = 1, alpha = 0.85) +
  pal_rank_color +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "#B40426", linewidth = 0.7) +
  annotate("text",
           x = max(df_all_summary$rank) * 0.98, y = -0.52,
           label = "dep = −0.5",
           color = "#B40426", hjust = 1, vjust = 1.4,
           size = 4, fontface = "italic") +
  scale_x_continuous(
    breaks = c(1, seq(100, nrow(df_all_summary), by = 100)),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Gene rank",
    y = "Mean dependency score across ovarian cancer cell lines"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text  = element_text(size = 7.5),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 12)
  ) +
  guides(fill = "none")

# Zoom plot
p_zoom <- ggplot(zoom_data, aes(rank, mean_dep, color = mean_dep)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean_dep - sd_dep, ymax = mean_dep + sd_dep),
                width = 0.2, linewidth = 0.4, alpha = 0.7) +
  geom_label_repel(
    aes(label = gene),
    size = 3.5, fontface = "bold",
    fill = "white", color = "black",
    box.padding = 0.3, point.padding = 0.2,
    segment.color = "#B40426",
    segment.size = 0.3,
    max.overlaps = 20, seed = 42
  ) +
  scale_color_gradient2(
    low = "#B40426", mid = "#FFDD57", high = "#2166AC",
    midpoint = 0, guide = "none"
  ) +
  scale_x_continuous(breaks = 1:TOP_N) +
  labs(
    x = NULL,
    y = NULL,
    title = "Top 10 most essential genes"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5),
    axis.title = element_blank(),   # 🔥 opcional extra (doble seguro)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = "grey70")
  )

# Combine plots
p_combined <- p_main +
  inset_element(
    p_zoom,
    left = 0.50, bottom = 0.09,
    right = 0.98, top = 0.56,
    align_to = "plot"
  )
p_combined
# ============================================================
# Save combined plot (correct way)
# ============================================================

Cairo::CairoPDF(
  file = paste0(OUT_DIR, "6_1_2_RankedProfile_allgenes.pdf"),
  width = 8, height = 8.5
)
print(p_combined)
dev.off()


png(paste0(OUT_DIR, "6_1_2_RankedProfile_allgenes.png"),
    width = 8, height = 8.5, units = "in", res = 600)
print(p_combined)
dev.off()


#save.image("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/5_3_HM_DepMap.RData")

