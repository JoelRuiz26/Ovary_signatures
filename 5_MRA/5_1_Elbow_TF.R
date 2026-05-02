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
  "~/Ovary_signatures/3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv") %>% filter(n_sources==3)
# TF list
TF_list <- vroom(
  "~/Ovary_signatures/5_MRA/core_mra.tsv",show_col_types=FALSE)

anot_genes <- readRDS("~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_2_gene_annotation_UNION_universe.rds")

# core genes
core_list <- full_genes %>%
  pull(gene)
length(core_list) #[1] 648

# TF MRA
df <- TF_list %>%
  filter(padj_emp < 0.05) %>%
  mutate(abs_NES = abs(NES)) %>%
  arrange(rank_abs_NES)
df_NES2 <- df %>% filter(abs_NES >= 2)

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

elbow_idx <- which.max(distances) #[1] 68
elbow_point <- df[elbow_idx, ]

p_elbow <- ggplot(df, aes(rank_abs_NES, abs_NES)) +
  geom_line(color = "black", linewidth = 0.7) +
  
  # zona seleccionada (core regulons)
  annotate("rect",
           xmin = 0, xmax = elbow_idx,
           ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.08) +
  
  # punto del elbow
  geom_point(data = elbow_point,
             size = 2.5, color = "#D55E00") +
  
  # etiqueta discreta
  geom_text(data = elbow_point,
            aes(label = elbow_idx),
            vjust = -1, size = 3.5, color = "black") +
  
  labs(x = "Rank", y = "|NES|") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_elbow


# Gene Name symbol
tf_names <- df %>%
  filter(rank_abs_NES <= elbow_point$rank_abs_NES) %>%
  mutate(TF = toupper(TF)) %>% 
  dplyr::select(TF, NES, direction) %>% as.data.frame()


saveRDS(tf_names,"~/Ovary_signatures/5_MRA/5_0_1_Top_elbow_regulones.rds")


tf_names <- tf_names%>%
  pull(TF) %>%
  unique()

length(tf_names) #328


#Gene enbl
#tf_list <- anot_genes %>%
#  filter(toupper(gene) %in% tf_names) %>%
#  pull(ensembl)
#length(tf_list) #[1] 320

# ------------------------------------------------------------
# Plot elbow
# ------------------------------------------------------------

outdir <- file.path("~/Ovary_signatures/5_MRA/")


ggsave(file.path(outdir, "5_0_Elbow_plot_regulons.pdf"),
       p_elbow, width = 6, height = 4)

ggsave(file.path(outdir, "5_0_Elbow_plot_regulons.png"),
       p_elbow, width = 6, height = 4, dpi = 600)



# ------------------------------------------------------------
# b) Core genes in regulons
# ------------------------------------------------------------

Master_Regulator <- readRDS("~/Ovary_signatures/5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds")

# Filtrar regulones
regulons_filtered <- Master_Regulator[names(Master_Regulator) %in% tf_names]
saveRDS(regulons_filtered,"~/Ovary_signatures/5_MRA/5_0_Top_elbow_regulones.rds")


# Extraer todos los genes de regulons_filtered
regulon_genes <- unique(unlist(lapply(regulons_filtered, function(x) {
  names(x$tfmode)
})))

# Guardar pool de genes
saveRDS(regulon_genes,
        "~/Ovary_signatures/5_MRA/5_1_regulon_genes_pool.rds")


# Filtrar core_genes con ese pool
core_genes_filtered <- full_genes %>%
  filter(toupper(gene) %in% toupper(regulon_genes))

# Guardar core genes filtrados
saveRDS(core_genes_filtered,
        "~/Ovary_signatures/5_MRA/5_2_core_genes_filtered.rds")


