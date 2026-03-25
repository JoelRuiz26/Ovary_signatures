#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(depmap)
  library(ggplot2)
  library(effsize)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 1) Load DGE
# ============================================================

DGE_PATH <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv"

dge <- vroom(DGE_PATH, show_col_types = FALSE) %>%
  filter(n_sources == 3,
         !is.na(mean_log2FC)) %>%
  mutate(gene = toupper(gene)) %>%
  group_by(gene) %>%
  summarise(median_log2FC = median(mean_log2FC, na.rm = TRUE), .groups="drop")

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

# ============================================================
# 6) Statistics
# ============================================================

wilcox_res <- wilcox.test(median_dep ~ direction, data = df, exact = FALSE)

cliff_res <- cliff.delta(
  df$median_dep[df$direction == "Upregulated"],
  df$median_dep[df$direction == "Downregulated"]
)

pearson_res <- cor.test(df$median_log2FC, df$median_dep, method = "pearson")

# ============================================================
# 7) Plot
# ============================================================
p <- ggplot(df, aes(x = median_log2FC, y = median_dep)) +
  
  # puntos (monocromáticos)
  geom_point(color = "grey30", alpha = 0.7, size = 2) +
  
  # líneas de cuadrantes
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  
  # regresión (Pearson)
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  
  # tema limpio tipo publicación
  theme_classic(base_size = 11) +
  
  labs(
    title = "Association between gene expression changes and CRISPR dependency",
    subtitle = paste0(
      "Wilcoxon p < 2e⁻¹⁶ | Cliff’s δ ", round(cliff_res$estimate, 3),
      "\nPearson p ",format.pval(pearson_res$p.value, digits = 2), 
      " | r ", round(pearson_res$estimate, 3)
    ),
    x = "log2 Fold Change",
    y = "CRISPR dependency score in cell lines"
  )

print(p)
# ============================================================
# 8) Save
# ============================================================

ggsave("~/Ovary_signatures/6_Depmap_ovary/quadrant_plot_publication.png",
       p, width = 6, height = 4, dpi = 600)



