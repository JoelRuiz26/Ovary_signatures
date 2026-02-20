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
# 2) Compute per-gene "noise terms"
#    Under DESeq2 NB model: Var = mu + alpha * mu^2
#    - Poisson part (mu) is the closest proxy to "technical noise"
#    - frac_technical = mu / (mu + alpha * mu^2)
# =========================
compute_noise <- function(dds) {
  mu   <- rowMeans(counts(dds, normalized = TRUE))
  disp <- dispersions(dds)  # alpha
  
  # split variance into components
  var_tech <- mu
  var_bio  <- disp * mu^2
  var_tot  <- var_tech + var_bio
  
  tibble(
    gene = rownames(dds),
    mean_expression = mu,
    dispersion = disp,
    bcv = sqrt(disp),
    tech_noise = sqrt(var_tech),                 # "expected sampling noise" (SD scale)
    frac_technical = var_tech / var_tot          # fraction of variance explained by Poisson term
  ) %>%
    filter(is.finite(mean_expression), mean_expression > 0,
           is.finite(dispersion), dispersion > 0,
           is.finite(frac_technical), frac_technical >= 0, frac_technical <= 1) %>%
    arrange(desc(frac_technical))
}

noise_GTEX <- compute_noise(dds_tumor_gtex)
noise_AE   <- compute_noise(dds_tumor_ae)

# =========================
# "Technical noise" comparison (paired Wilcoxon)
# using frac_technical per gene (aligned by gene)
# =========================
df_tech_cmp <- inner_join(
  noise_AE   %>% select(gene, frac_technical),
  noise_GTEX %>% select(gene, frac_technical),
  by = "gene",
  suffix = c("_AE", "_GTEx")
)

wilcox_tech <- wilcox.test(
  df_tech_cmp$frac_technical_AE,
  df_tech_cmp$frac_technical_GTEx,
  paired = TRUE
)

summary_tech <- df_tech_cmp %>%
  summarise(
    median_AE = median(frac_technical_AE),
    median_GTEx = median(frac_technical_GTEx),
    median_ratio = median(frac_technical_AE / frac_technical_GTEx)
  )

wilcox_tech
summary_tech

# =========================
# 3) Map gene subsets (same as you did)
# =========================
all_genes <- vroom::vroom("~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv")

map_src <- all_genes %>%
  dplyr::select(ensembl, sources) %>%
  mutate(dataset = toupper(gsub(",", "_", sources))) %>%
  distinct(ensembl, dataset)

noise_GTEX <- noise_GTEX %>% left_join(map_src, by = c("gene" = "ensembl"))
noise_AE   <- noise_AE   %>% left_join(map_src, by = c("gene" = "ensembl"))

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
# 1) GLOBAL plot (GTEx vs AE)  [Technical noise share]
# -------------------------
df_plot_global <- bind_rows(
  noise_GTEX %>% mutate(dataset = "GTEx"),
  noise_AE   %>% mutate(dataset = "AE")
) %>%
  filter(is.finite(frac_technical), frac_technical >= 0, frac_technical <= 1) %>%
  filter(dataset != "GEO")

dx_global <- (max(df_plot_global$frac_technical, na.rm = TRUE) -
                min(df_plot_global$frac_technical, na.rm = TRUE)) / 120

p_global <- ggplot(df_plot_global, aes(x = frac_technical, fill = dataset, color = dataset)) +
  geom_density(
    aes(y = after_stat(density) * 100 * dx_global),
    alpha = 0.30,
    linewidth = 1.0,
    adjust = 1.2
  ) +
  scale_color_manual(values = col_global) +
  scale_fill_manual(values = col_global) +
  labs(
    title = "Fraction of model variance attributable to Poisson sampling",
    x = "Poisson sampling fraction of variance",
    y = "Genes (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

print(p_global + coord_cartesian(xlim = c(0, 0.03)))


# -------------------------
# 2) SUBSET plot (LINES ONLY)  [Technical noise share]
# -------------------------
df_plot_subset <- bind_rows(
  noise_GTEX %>% mutate(list = "GTEx"),
  noise_AE   %>% mutate(list = "AE")
) %>%
  filter(!is.na(dataset), is.finite(frac_technical), frac_technical >= 0, frac_technical <= 1) %>%
  filter(dataset != "GEO") %>%
  filter(
    (list == "GTEx" & grepl("GTEX", dataset)) |
      (list == "AE"  & grepl("AE",   dataset))
  )

dx_subset <- (max(df_plot_subset$frac_technical, na.rm = TRUE) -
                min(df_plot_subset$frac_technical, na.rm = TRUE)) / 120

p_subset <- ggplot(df_plot_subset, aes(x = frac_technical)) +
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
    breaks = c("AE","AE_GEO","GTEX","GTEX_GEO","GTEX_AE","GTEX_AE_GEO")
  ) +
  scale_linetype_manual(values = c("GTEx" = "dashed", "AE" = "solid")) +
  labs(
    title = "Fraction of model variance attributable to Poisson sampling",
    x = "Poisson sampling fraction of variance",
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

print(p_subset + coord_cartesian(xlim = c(0, 0.01)))


# -------------------------
# 3) SUBSET plot (COMBINED AE + GTEx)  [Technical noise share]
# -------------------------
df_plot_subset_combined <- bind_rows(
  noise_GTEX,
  noise_AE
) %>%
  filter(!is.na(dataset), is.finite(frac_technical), frac_technical >= 0, frac_technical <= 1) %>%
  filter(dataset != "GEO")

dx_subset_combined <- (max(df_plot_subset_combined$frac_technical, na.rm = TRUE) -
                         min(df_plot_subset_combined$frac_technical, na.rm = TRUE)) / 120

p_subset_combined <- ggplot(df_plot_subset_combined, aes(x = frac_technical)) +
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
    breaks = c("AE","AE_GEO","GTEX","GTEX_GEO","GTEX_AE","GTEX_AE_GEO")
  ) +
  labs(
    title = "Fraction of model variance attributable to Poisson sampling",
    x = "Poisson sampling fraction of variance",
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

print(p_subset_combined + coord_cartesian(xlim = c(0, 0.01)))

#####

full_df <- bind_rows(
  noise_GTEX %>% mutate(source = "GTEx"),
  noise_AE   %>% mutate(source = "AE")
) %>%
  filter(!is.na(dataset), dataset != "GEO")
kruskal.test(frac_technical ~ dataset, data = full_df)

