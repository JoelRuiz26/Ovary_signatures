#load("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/6_2_Image_Ovary_correlation_Depmap_DEG.RData")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(ggplot2)
  library(effsize)
  library(ggplot2)
  library(ggpubr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 1) Load DGE
# ============================================================

DGE_PATH <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv"

dge <- vroom(DGE_PATH, show_col_types = FALSE) %>%
  filter(n_sources == 3,
         !is.na(log2FC_ae),
         abs(log2FC_ae) > 1) %>%
  mutate(gene = toupper(gene)) %>%
  group_by(gene) %>%
  summarise(median_log2FC = median(log2FC_ae, na.rm = TRUE), .groups="drop")
#  summarise(median_log2FC = median(mean_log2FC, na.rm = TRUE), .groups="drop")

# ============================================================
# 2) DepMap ovarian
# ============================================================

crispr <- depmap_crispr() %>%
  filter(grepl("_OVARY$", cell_line)) %>%
  mutate(
    gene = toupper(gene_name),
    dep = dependency
  )

# ============================================================
# 3) Median dependency
# ============================================================

dep_med <- crispr %>%
  group_by(gene) %>%
  summarise(median_dep = median(dep, na.rm = TRUE), .groups="drop")
summary(dep_med$median_dep)
# ============================================================
# 4) Merge
# ============================================================

df <- inner_join(dge, dep_med, by = "gene")

# ============================================================
# 5) Define groups
# ============================================================

df <- df %>%
  mutate(
    direction = ifelse(median_log2FC > 0, "Upregulated", "Downregulated")
  )
vroom_write(df, "~/Ovary_signatures/6_Depmap_ovary/6_2_0_Dependency_DEG_core.tsv")


# ============================================================
# 6) Statistics
# ============================================================

wilcox_res <- wilcox.test(median_dep ~ direction, data = df, exact = FALSE)

cliff_res <- cliff.delta(
  df$median_dep[df$direction == "Upregulated"],
  df$median_dep[df$direction == "Downregulated"]
)

pearson_res <- cor.test(df$median_log2FC, df$median_dep, method = "pearson")
spearman_res <- cor.test(df$median_log2FC, df$median_dep, method = "spearman")
# ============================================================
# 7) Plot
# ============================================================
p <- ggplot(df, aes(x = median_log2FC, y = median_dep)) +
  
  # puntos (monocromáticos)
  geom_point(color = "grey30", alpha = 0.7, size = 2) +
  
  # líneas de cuadrantes
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  
  # regresión
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  
  # tema limpio tipo publicación
  theme_classic(base_size = 12) +
  
  labs(
    title = "Association between gene expression changes and CRISPR dependency",
    subtitle = paste0(
      "\nSpearman r ", round(spearman_res$estimate, 3),
      " | p ", format.pval(spearman_res$p.value, digits = 3)
    ),
    x = "log2 Fold Change",
    y = "Median CRISPR dependency score in cell lines"
  )

print(p)

#####Wilcoxon###
# ============================================================
#Wilcoxon plot
# ============================================================

# Reordenar factor para estética
df$direction <- factor(df$direction, levels = c("Downregulated", "Upregulated"))

p_wilcox <- ggplot(df, aes(x = direction, y = median_dep, fill = direction)) +
  
  geom_violin(trim = FALSE, scale = "width", alpha = 0.6, color = NA) +
  
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black", linewidth = 0.4) +
  
  geom_point(position = position_jitter(width = 0.1), alpha = 0.2, size = 0.7) +
  
  scale_fill_manual(values = c("#3C5488FF", "#B09C85FF")) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  
  coord_cartesian(ylim = c(-2.5, 0.5)) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 12),
    
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12),
    
    # 🔥 ESTO ES LO QUE QUIERES
    axis.text.x = element_text(face = "bold", size = 13)
  ) +
  labs(
    title = NULL,
    subtitle = paste0(
      "Wilcoxon p = ", format.pval(wilcox_res$p.value, digits = 3),
      " | Cliff’s δ = ", round(cliff_res$estimate, 2)
    ),
    y = "Mean dependency score across ovarian cancer cell lines"
  )

print(p_wilcox)

# ============================================================
# 8) Save
# ============================================================

ggsave("~/Ovary_signatures/6_Depmap_ovary/6_2_1_FourCUadtrant_plot.png",
       p, width = 5, height = 4.5, dpi = 600)

# ============================================================
# Save Wilcoxon plot
# ============================================================

ggsave("~/Ovary_signatures/6_Depmap_ovary/6_2_2_wilcoxon_violin.png",
  p_wilcox, width = 5, height = 6.5, dpi = 600)
#save.image("/STORAGE/csbig/jruiz/Ovary_data/5_Depmap_ovary/6_2_Image_Ovary_correlation_Depmap_DEG.RData")
