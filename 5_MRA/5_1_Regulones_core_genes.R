### Identify regulons candidates for survival analysis
### a) Filter by elbow method 
### b) Identify the proportion of core genes present in Regulons (get ranket list of regulons)
### C) Get list of regulon(es) with higuest NES and enriched with core genes

library(dplyr)
library(tibble)
library(purrr)
library(vroom)
library(ggplot2)

# ------------------------------------------------------------
#  data
# ------------------------------------------------------------

full_genes <- vroom(
  "~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv") %>% filter(n_sources==3)
# TF list
TF_list <- vroom(
  "~/Ovary_signatures/5_MRA/core_mra.tsv",show_col_types=FALSE)

anot_genes <- readRDS("~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_2_gene_annotation_UNION_universe.rds")

# core genes
core_list <- full_genes %>%
  pull(ensembl)
length(core_list) #[1] 648

# TF MRA
df <- TF_list %>%
  filter(padj_emp < 0.05) %>%
  mutate(abs_NES = abs(NES)) %>%
  arrange(rank_abs_NES)

# ------------------------------------------------------------
# a) Filter by elbow method 
# ------------------------------------------------------------
# Elbow method
distances <- abs((df$abs_NES[nrow(df)] - df$abs_NES[1]) * df$rank_abs_NES -
                   (df$rank_abs_NES[nrow(df)] - df$rank_abs_NES[1]) * df$abs_NES +
                   df$rank_abs_NES[nrow(df)] * df$abs_NES[1] -
                   df$abs_NES[nrow(df)] * df$rank_abs_NES[1]) /
  sqrt((df$abs_NES[nrow(df)] - df$abs_NES[1])^2 +
         (df$rank_abs_NES[nrow(df)] - df$rank_abs_NES[1])^2)

elbow_idx <- which.max(distances)
elbow_point <- df[elbow_idx, ]

p_elbow <- ggplot(df, aes(rank_abs_NES, abs_NES)) +
  geom_line() +
  geom_vline(xintercept = elbow_point$rank_abs_NES,
             linetype = "dashed", color = "red") +
  geom_text(data = elbow_point,
            aes(label = rank_abs_NES),
            vjust = -1, color = "black") +
  theme_minimal()


# Gene Name symbol
tf_names <- df %>%
  filter(rank_abs_NES <= elbow_point$rank_abs_NES) %>%
  mutate(TF = toupper(TF)) %>%
  pull(TF) %>%
  unique()

length(tf_names) #[1] 328

#Gene enbl
#tf_list <- anot_genes %>%
#  filter(toupper(gene) %in% tf_names) %>%
#  pull(ensembl)
#length(tf_list) #[1] 320

# ------------------------------------------------------------
# b) Core genes in regulons
# ------------------------------------------------------------

Master_Regulator <- readRDS("/STORAGE/csbig/jruiz/Ovary_data/6_MRA/cancer_ovary_regulon_300bt_p1e-8.rds")

# Filtrar regulones
regulons_filtered <- Master_Regulator[names(Master_Regulator) %in% tf_names]
saveRDS(regulons_filtered,"/STORAGE/csbig/jruiz/Ovary_data/6_MRA/6_0_Top_elbow_regulones.rds")


# Métricas
regulon_stats <- map_dfr(names(regulons_filtered), function(tf){
  genes <- names(regulons_filtered[[tf]]$tfmode)
  tibble(
    TF = tf,
    n_genes = length(genes),
    n_core = sum(genes %in% full_genes$gene),
    prop_core = n_core / n_genes
  )
})

# ------------------------------------------------------------
# TOP número
# ------------------------------------------------------------

top30_n <- regulon_stats %>%
  slice_max(n_core, n = 32, with_ties = FALSE)

top5_n <- regulon_stats %>%
  slice_max(n_core, n = 32, with_ties = FALSE)

# plot
p_n <- ggplot(top30_n, aes(x = reorder(TF, n_core))) +
  geom_col(aes(y = n_genes), fill = "grey80") +
  geom_col(aes(y = n_core), fill = "steelblue") +
  geom_text(aes(y = n_core, label = n_core),
            hjust = -0.2, size = 3) +
  coord_flip() +
  theme_minimal()

# guardar regulones top5
top5_n_regulons <- map(
  regulons_filtered[names(regulons_filtered) %in% top5_n$TF],
  function(reg){
    keep <- names(reg$tfmode) %in% full_genes$gene
    list(
      tfmode = reg$tfmode[keep],
      likelihood = reg$likelihood[keep]
    )
  }
)
# ------------------------------------------------------------
# TOP proporción
# ------------------------------------------------------------

top30_prop <- regulon_stats %>%
  arrange(desc(prop_core), desc(n_core)) %>%
  slice_head(n = 32)

top5_prop <- regulon_stats %>%
  arrange(desc(prop_core), desc(n_core)) %>%
  slice_head(n = 32)

# plot
p_prop <- ggplot(top30_prop, aes(x = reorder(TF, prop_core))) +
  geom_col(aes(y = 1), fill = "grey80") +
  geom_col(aes(y = prop_core), fill = "darkgreen") +
  geom_text(aes(
    y = prop_core,
    label = paste0(round(prop_core * 100, 1), "% (", n_core, ")")
  ),
  hjust = -0.2, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()


# guardar regulones top5
top5_prop_regulons <- map(
  regulons_filtered[names(regulons_filtered) %in% top5_prop$TF],
  function(reg){
    keep <- names(reg$tfmode) %in% full_genes$gene
    list(
      tfmode = reg$tfmode[keep],
      likelihood = reg$likelihood[keep]
    )
  }
)

# ------------------------------------------------------------
# Plot número
# ------------------------------------------------------------
outdir <- file.path("/STORAGE/csbig/jruiz/Ovary_data/6_MRA/")
ggsave(file.path(outdir, "6_1_Top30_n_core.pdf"), p_n, width = 7, height = 10)

ggsave(file.path(outdir, "6_1_Top30_n_core.png"), p_n,
       width = 7, height = 10, dpi = 600)

saveRDS(top5_n_regulons,
        file = file.path(outdir, "6_1_top5_regulons_n_core.rds"))

# ------------------------------------------------------------
# Plot proporción
# ------------------------------------------------------------

ggsave(file.path(outdir, "6_2_Top30_prop_core.pdf"), p_prop, width = 7, height = 10)

ggsave(file.path(outdir, "6_2_Top30_prop_core.png"), p_prop,
       width = 7, height = 10, dpi = 600)

saveRDS(top5_prop_regulons,
        file = file.path(outdir, "6_2_top5_regulons_prop_core.rds"))

# ------------------------------------------------------------
# Plot elbow
# ------------------------------------------------------------

ggsave(file.path(outdir, "6_0_Elbow_plot_regulons.pdf"),
       p_elbow, width = 6, height = 4)

ggsave(file.path(outdir, "6_0_Elbow_plot_regulons.png"),
       p_elbow, width = 6, height = 4, dpi = 600)


