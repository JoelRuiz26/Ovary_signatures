#load("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/5_3_HM_DepMap.RData")


# ============================================================
# 6) TF regulatory activity vs CRISPR dependency
# ============================================================

tf_dep_all <- crispr_ovary_all %>%
  group_by(gene_u) %>%
  summarise(
    mean_dependency = median(dep_scaled_m11, na.rm=TRUE),
    .groups="drop"
  ) %>%
  dplyr::rename(TF=gene_u)

tf_plot <- TF_signif %>%
  dplyr::transmute(TF=toupper(TF), NES, rank_abs_NES) %>%
  left_join(tf_dep_all, by="TF") %>%
  filter(!is.na(mean_dependency)) %>%
  mutate(
    quadrant = case_when(
      NES > 0 & mean_dependency < 0 ~ "Activated & essential",
      NES < 0 & mean_dependency < 0 ~ "Repressed but essential",
      NES > 0 & mean_dependency > 0 ~ "Activated non-essential",
      NES < 0 & mean_dependency > 0 ~ "Repressed non-essential"
    )
  )

top_label <- tf_plot %>%
  arrange(mean_dependency) %>%
  slice_head(n=20)

p_quad <- ggplot(tf_plot, aes(NES, mean_dependency)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") +
  geom_point(aes(color=quadrant), size=2.5, alpha=0.8) +
  geom_text_repel(data=top_label, aes(label=TF), size=3.5, max.overlaps=Inf) +
  scale_color_manual(values=c(
    "Activated & essential"="#D73027",
    "Repressed but essential"="#4575B4",
    "Activated non-essential"="#FDAE61",
    "Repressed non-essential"="#74ADD1"
  )) +
  labs(
    x="Regulatory activity (NES)",
    y="CRISPR dependency (mean scaled)",
    color="TF class",
    title="Transcription factor regulatory activity in tumors and dependency in ovarian cell line"
  ) +
  theme_bw(base_size=14)

print(p_quad)

pdf("~/Ovary_signatures/6_Depmap_ovary/6_1_2_TF_NES_vs_dependency.pdf",8,7)
print(p_quad); dev.off()

png("~/Ovary_signatures/6_Depmap_ovary/6_1_2_TF_NES_vs_dependency.png",
    width=8,height=7,units="in",res=600)
print(p_quad); dev.off()


####
# ============================================================
# 7) Heatmap TF (Top 40 |NES|) highlighting consensus genes
# ============================================================

# Top 40 TF by absolute NES
tf_top40 <- TF_signif %>%
  arrange(rank_abs_NES) %>%
  slice_head(n = 80) %>%
  pull(TF)

# Keep only TF present in DepMap
tf_top40 <- intersect(tf_top40, unique(crispr_ovary_all$gene_u))

# Build TF dependency matrix
df_tf <- crispr_ovary_all %>%
  filter(gene_u %in% tf_top40) %>%
  transmute(
    gene = gene_u,
    cell = sub("_OVARY$", "", cell_line),
    dep_scaled_m11
  ) %>%
  distinct()

mat_tf <- df_tf %>%
  pivot_wider(names_from = cell, values_from = dep_scaled_m11) %>%
  column_to_rownames("gene") %>%
  as.matrix()

mat_tf[is.na(mat_tf)] <- NA

# Order by essentiality
mat_tf <- mat_tf[order(rowMeans(mat_tf)), ]

# Highlight TF also present in consensus genes
cons_gene_names <- cons_genes$gene_u

is_cons <- rownames(mat_tf) %in% cons_gene_names

row_colors_tf <- ifelse(is_cons, "#228B22", "black")
row_faces_tf  <- ifelse(is_cons, "bold", "plain")

# Dendrogram
row_dend_tf <- as.dendrogram(hclust(dist(mat_tf)))
row_dend_tf <- reorder(row_dend_tf, wts = -rowMeans(mat_tf))
row_dend_tf <- rev(row_dend_tf)

# Heatmap
ht_tf_topNES <- Heatmap(
  mat_tf,
  name="Scaled\ndependency",
  col=col_fun,
  cluster_rows=row_dend_tf,
  cluster_columns=TRUE,
  show_row_names=TRUE,
  show_column_names=TRUE,
  row_names_gp=gpar(
    fontsize=8,
    col=row_colors_tf,
    fontface=row_faces_tf
  ),
  column_names_gp=gpar(fontsize=8),
  column_title="Transcription Factor Dependency Across Ovarian Cancer Cell Lines",
  column_title_gp=gpar(fontsize=20,fontface="bold"),
  rect_gp=gpar(col=NA),
  border=FALSE
)

grid.newpage()
draw(ht_tf_topNES, heatmap_legend_side="right")

# Save
pdf("~/Ovary_signatures/6_Depmap_ovary/6_1_3_Heatmap_TF_topNES.pdf",
    width=7.5, height=6.5)
grid.newpage()
draw(ht_tf_topNES, heatmap_legend_side="right")
dev.off()

png("~/Ovary_signatures/6_Depmap_ovary/6_1_3_Heatmap_TF_topNES.png",
    width=7.5, height=6.5,
    units="in", res=600)
grid.newpage()
draw(ht_tf_topNES, heatmap_legend_side="right")
dev.off()
