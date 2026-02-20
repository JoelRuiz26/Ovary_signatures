suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(ggplot2)
  library(DESeq2)
})

OUT_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor"

# =========================
# 1) Load DESeq2 objects (NB model already fit)
# =========================
dds_tumor_gtex <- readRDS(file.path(OUT_DIR, "3_3_2_dds_GTEX.rds"))
dds_tumor_ae   <- readRDS(file.path(OUT_DIR, "3_3_2_dds_AE.rds"))

# =========================
# 2) Compute NB variability per gene (DESeq2 dispersion)
#    dispersion = alpha in Var = mu + alpha * mu^2
# =========================
compute_nb_var <- function(dds) {
  mu <- rowMeans(counts(dds, normalized = TRUE))
  disp <- dispersions(dds)
  
  tibble(
    gene = rownames(dds),
    mean_expression = mu,
    dispersion = disp,
    bcv = sqrt(disp)
  ) %>%
    filter(is.finite(mean_expression), mean_expression > 0,
           is.finite(dispersion), dispersion > 0) %>%
    arrange(desc(bcv))
}

nb_GTEX <- compute_nb_var(dds_tumor_gtex)
nb_AE   <- compute_nb_var(dds_tumor_ae)


# =========================
# BCV comparison (paired Wilcoxon)
# =========================
df_bcv_cmp <- inner_join(
  nb_AE %>% dplyr::select(gene, bcv),
  nb_GTEX %>% dplyr::select(gene, bcv),
  by = "gene",
  suffix = c("_AE", "_GTEx")
)

wilcox_bcv <- wilcox.test(
  df_bcv_cmp$bcv_AE,
  df_bcv_cmp$bcv_GTEx,
  paired = TRUE
)

summary_stats <- df_bcv_cmp %>%
  summarise(
    median_AE = median(bcv_AE),
    median_GTEx = median(bcv_GTEx),
    median_ratio = median(bcv_AE / bcv_GTEx)
  )

wilcox_bcv
summary_stats

#median_AE median_GTEx median_ratio
#<dbl>       <dbl>        <dbl>
#  1     0.710       0.661         1.05

#AE tiene aproximadamente 5% mayor BCV mediano que GTEx


# =========================
# 3) Map gene subsets (same as you did)
# =========================
all_genes <- vroom::vroom("~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv")

map_src <- all_genes %>%
  dplyr::select(ensembl, sources) %>%
  mutate(dataset = toupper(gsub(",", "_", sources))) %>%
  distinct(ensembl, dataset)

nb_GTEX <- nb_GTEX %>% left_join(map_src, by = c("gene" = "ensembl"))
nb_AE   <- nb_AE   %>% left_join(map_src, by = c("gene" = "ensembl"))



# =========================
# COLORS (manual, matching your palette)
# =========================
col_global <- c(
  "GTEx" = "#7A8DB8",
  "AE"   = "#F4D96B"
)

col_subset <- c(
  "AE"          = "#F4D96B",  
  "GTEX"        = "#7A8DB8",  
  "GTEX_AE"     = "#5AAE61",  
  "AE_GEO"      = "#F9EDB2",  
  "GTEX_GEO"    = "#B5C3E3",  
  "GTEX_AE_GEO" = "#E41A1C" 
)

# =========================
# PLOTS
# =========================

# -------------------------
# 1) GLOBAL plot (GTEx vs AE)  [BCV]
# -------------------------
df_plot_global <- bind_rows(
  nb_GTEX %>% mutate(dataset = "GTEx"),
  nb_AE   %>% mutate(dataset = "AE")
) %>%
  filter(is.finite(bcv), bcv > 0) %>%
  filter(dataset != "GEO")

dx_global <- (max(df_plot_global$bcv, na.rm = TRUE) -
                min(df_plot_global$bcv, na.rm = TRUE)) / 120

p_global <- ggplot(df_plot_global, aes(x = bcv, fill = dataset, color = dataset)) +
  geom_density(
    aes(y = after_stat(density) * 100 * dx_global),
    alpha = 0.30,
    linewidth = 1.0,
    adjust = 1.2
  ) +
  scale_color_manual(values = col_global) +
  scale_fill_manual(values = col_global) +
  labs(
    title = "Biological Coefficient Variation distribution",
    x = "BCV = sqrt(dispersion)",
    y = "Genes (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

print(p_global + coord_cartesian(xlim = c(0, 2)))


# -------------------------
# 2) SUBSET plot (LINES ONLY)  [BCV]
# -------------------------
df_plot_subset <- bind_rows(
  nb_GTEX %>% mutate(list = "GTEx"),
  nb_AE   %>% mutate(list = "AE")
) %>%
  filter(!is.na(dataset), is.finite(bcv), bcv > 0) %>%
  filter(dataset != "GEO") %>%
  filter(
    (list == "GTEx" & grepl("GTEX", dataset)) |
      (list == "AE"  & grepl("AE",   dataset))
  )

dx_subset <- (max(df_plot_subset$bcv, na.rm = TRUE) -
                min(df_plot_subset$bcv, na.rm = TRUE)) / 120

p_subset <- ggplot(df_plot_subset, aes(x = bcv)) +
  geom_density(
    aes(
      y = after_stat(density) * 100 * dx_subset,
      color = dataset,
      linetype = list,
      group = interaction(dataset, list)
    ),
    linewidth = 1.2,
    adjust = 1.2
  ) +
  scale_color_manual(
    values = col_subset,
    breaks = c(
      "AE",
      "AE_GEO",
      "GTEX",
      "GTEX_GEO",
      "GTEX_AE",
      "GTEX_AE_GEO"
    )
  )+
  scale_linetype_manual(values = c("GTEx" = "dashed", "AE" = "solid")) +
  labs(
    title = "Biological Coefficient Variation distributions by gene subset",
    x = "BCV = sqrt(dispersion)",
    y = "Genes (%)",
    color = "Gene subset",
    linetype = "List"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

print(p_subset + coord_cartesian(xlim = c(0, 2)))


#######################
# -------------------------
# 3) SUBSET plot (COMBINED AE + GTEx)  [BCV]
# -------------------------

df_plot_subset_combined <- bind_rows(
  nb_GTEX,
  nb_AE
) %>%
  filter(!is.na(dataset), is.finite(bcv), bcv > 0) %>%
  filter(dataset != "GEO")

dx_subset_combined <- (max(df_plot_subset_combined$bcv, na.rm = TRUE) -
                         min(df_plot_subset_combined$bcv, na.rm = TRUE)) / 120

p_subset_combined <- ggplot(df_plot_subset_combined, aes(x = bcv)) +
  geom_density(
    aes(
      y = after_stat(density) * 100 * dx_subset_combined,
      color = dataset,
      group = dataset
    ),
    linewidth = 1.3,
    adjust = 1.2
  ) +
  scale_color_manual(
    values = col_subset,
    breaks = c(
      "AE",
      "AE_GEO",
      "GTEX",
      "GTEX_GEO",
      "GTEX_AE",
      "GTEX_AE_GEO"
    )
  )+
  labs(
    title = "Biological Coefficient Variation distributions by gene subset",
    x = "BCV = sqrt(dispersion)",
    y = "Genes (%)",
    color = "Gene subset"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

print(p_subset_combined + coord_cartesian(xlim = c(0, 2)))
