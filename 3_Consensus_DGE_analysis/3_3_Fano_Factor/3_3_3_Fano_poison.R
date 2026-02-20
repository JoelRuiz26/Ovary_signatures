#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(ggplot2)
})

OUT_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor"

# =========================
# 1) Load normalized matrices
# =========================
norm_all_ovary <- readRDS(file.path(OUT_DIR, "3_3_1_expr_raw_Tumor_Ctlhomolog.rds"))
norm_all_auto  <- readRDS(file.path(OUT_DIR, "3_3_1_expr_raw_Tumor_CtlAuto.rds"))

# =========================
# 2) Function to compute Fano factor per gene
#    F = var(x)/mean(x)
#    computed on RAW scale (no log)
# =========================
compute_fano <- function(mat) {
  gene_mean <- rowMeans(mat)
  gene_var  <- apply(mat, 1, var)
  
  fano <- gene_var / gene_mean
  
  tibble(
    gene = rownames(mat),
    mean_expression = gene_mean,
    variance = gene_var,
    fano_factor = fano
  ) %>%
    filter(is.finite(mean_expression), mean_expression > 0) %>%   # avoid division artifacts
    filter(is.finite(fano_factor), fano_factor >= 0) %>%          # keep valid
    arrange(desc(fano_factor))
}

# =========================
# 3) Compute Fano for both datasets
# =========================
fano_GTEX <- compute_fano(norm_all_ovary)
fano_AE   <- compute_fano(norm_all_auto)

# =========================
# 4) Save results
# =========================
saveRDS(fano_GTEX, file.path(OUT_DIR, "3_3_3_fano_GTEX_raw.rds"))
saveRDS(fano_AE,   file.path(OUT_DIR, "3_3_3_fano_AE_raw.rds"))

cat("\nSummary fano_AE$fano_factor:\n")
print(summary(fano_AE$fano_factor))
cat("\nSummary fano_GTEX$fano_factor:\n")
print(summary(fano_GTEX$fano_factor))

# =========================
# 5) Identify genes shared across datasets (gene subsets)
# =========================
all_genes <- vroom::vroom("~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_2_allgenes.tsv")

map_src <- all_genes %>%
  dplyr::select(ensembl, sources) %>%
  mutate(dataset = toupper(gsub(",", "_", sources))) %>%  # "GTEx,AE,GEO" -> "GTEX_AE_GEO"
  distinct(ensembl, dataset)

fano_GTEX <- fano_GTEX %>% left_join(map_src, by = c("gene" = "ensembl"))
fano_AE   <- fano_AE   %>% left_join(map_src, by = c("gene" = "ensembl"))

# =========================
# 6) Palettes
# =========================
palette_global <- c("GTEx" = "#7A8DB8", "AE" = "#F4D96B")

palette_dataset <- c(
  "AE" = "#F4D96B",
  "GTEX" = "#7A8DB8",
  "GTEX_AE" = "#5AAE61",
  "AE_GEO" = "#F9EDB2",
  "GTEX_GEO" = "#B5C3E3",
  "GTEX_AE_GEO" = "#E41A1C"
)

# =========================
# 7) PLOTS (LINEAR axis)
# =========================

# -------------------------
# 7.1) GLOBAL plot (GTEx vs AE)  [LINEAR]
# -------------------------
df_plot_global <- bind_rows(
  fano_GTEX %>% mutate(dataset_global = "GTEx"),
  fano_AE   %>% mutate(dataset_global = "AE")
) %>%
  filter(is.finite(fano_factor), fano_factor >= 0)

# Choose a linear zoom limit (keeps plot interpretable)
# 99.5% works well when there are extreme tails (millions)
xlim_global <- as.numeric(quantile(df_plot_global$fano_factor, probs = 0.995, na.rm = TRUE))

# IMPORTANT: compute dx on the VISIBLE range, not the full max (millions)
df_global_zoom <- df_plot_global %>%
  filter(fano_factor >= 0, fano_factor <= xlim_global)

dx_global <- (max(df_global_zoom$fano_factor, na.rm = TRUE) -
                min(df_global_zoom$fano_factor, na.rm = TRUE)) / 120

p_global <- ggplot(df_global_zoom, aes(x = fano_factor, fill = dataset_global, color = dataset_global)) +
  geom_density(
    aes(y = after_stat(density) * 100 * dx_global),
    alpha = 0.30,
    linewidth = 1.0,
    adjust = 1.2
  ) +
  scale_fill_manual(values = palette_global) +
  scale_color_manual(values = palette_global) +
  labs(
    title = "Fano factor distribution",
    x = "Fano factor",
    y = "Genes (%)"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# still use coord_cartesian to enforce the same viewing window
print(p_global + coord_cartesian(xlim = c(0, 500)))


# -------------------------
# 7.2) SUBSET plot (LINES ONLY)  [LINEAR]
# -------------------------
df_plot_subset <- bind_rows(
  fano_GTEX %>% mutate(list = "GTEx"),
  fano_AE   %>% mutate(list = "AE")
) %>%
  filter(!is.na(dataset), is.finite(fano_factor), fano_factor >= 0) %>%
  filter(
    (list == "GTEx" & grepl("GTEX", dataset)) |
      (list == "AE"  & grepl("AE",   dataset))
  )

# Slightly tighter zoom for subset plot (fine-grained differences)
xlim_subset <- as.numeric(quantile(df_plot_subset$fano_factor, probs = 0.99, na.rm = TRUE))

df_subset_zoom <- df_plot_subset %>%
  filter(fano_factor >= 0, fano_factor <= xlim_subset)

dx_subset <- (max(df_subset_zoom$fano_factor, na.rm = TRUE) -
                min(df_subset_zoom$fano_factor, na.rm = TRUE)) / 120

p_subset <- ggplot(df_subset_zoom, aes(x = fano_factor)) +
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
  scale_color_manual(values = palette_dataset, na.value = "grey40") +
  scale_linetype_manual(values = c("GTEx" = "dashed", "AE" = "solid")) +
  labs(
    title = "Fano factor distributions by gene subset",
    x = "Fano factor",
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

print(p_subset + coord_cartesian(xlim = c(0, 4000)))


#####Violin

df_violin <- df_subset_zoom %>%
  mutate(
    dataset = factor(
      dataset,
      levels = c("AE","GTEX","AE_GEO","GTEX_GEO","GTEX_AE","GTEX_AE_GEO")
    ),
    list = factor(list, levels = c("GTEx", "AE"))
  )

# --- plot base (NO cambies estética) ---
p_subset_violin <- ggplot(
  df_violin,
  aes(
    x = dataset,
    y = fano_factor,
    fill = dataset,
    group = interaction(dataset, list)
  )
) +
  geom_violin(
    trim      = TRUE,
    scale     = "width",
    linewidth = 0.3,
    position  = position_dodge(width = 0.85)
  ) +
  geom_boxplot(
    width         = 0.12,
    outlier.shape = NA,
    linewidth     = 0.25,
    fill          = NA,
    position      = position_dodge(width = 0.85)
  ) +
  scale_fill_manual(values = palette_dataset, na.value = "grey80") +
  labs(
    title = "Fano factor by gene subset",
    x = "Gene subset",
    y = "Fano factor"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# =========================
# Brackets + Δ
# =========================

# helper: x base (1..n) por nivel del factor dataset
x_base <- setNames(seq_along(levels(df_violin$dataset)), levels(df_violin$dataset))

# ----- A) Δ dentro de columna (dos violines) para GTEX_AE y GTEX_AE_GEO -----
med_intra <- df_violin %>%
  filter(dataset %in% c("GTEX_AE", "GTEX_AE_GEO")) %>%
  group_by(dataset, list) %>%
  summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
  mutate(x0 = x_base[as.character(dataset)])

# offset consistente con dodge(width=0.85): posición aprox de cada violín respecto al centro
offset <- 0.85 / 4

pair_intra <- med_intra %>%
  mutate(x = ifelse(list == "GTEx", x0 - offset, x0 + offset)) %>%
  group_by(dataset, x0) %>%
  summarise(
    y_low  = min(med),
    y_high = max(med),
    delta  = y_high - y_low,
    .groups = "drop"
  ) %>%
  mutate(
    col = "black",
    x_br = x0,                      # bracket centrado en la columna
    tick = 0.10,                    # tamaño de patitas
    y_text = y_high + 70,
    label = paste0("\u0394=", round(delta, 1))
  )

# ----- B) Δ entre columnas (AE_GEO vs GTEX_GEO) -----
med_inter <- df_violin %>%
  filter(dataset %in% c("AE_GEO", "GTEX_GEO")) %>%
  group_by(dataset) %>%
  summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
  mutate(x = x_base[as.character(dataset)])

# bracket vertical en el punto medio entre las dos columnas
pair_inter <- tibble(
  x_br  = mean(med_inter$x),
  y_low = min(med_inter$med),
  y_high= max(med_inter$med),
  delta = max(med_inter$med) - min(med_inter$med),
  col   = "black",
  tick  = 0.12,
  y_text = max(med_inter$med) + 70,
  label = paste0("\u0394=", round(delta, 1))
)

# --- Añadir brackets + texto Δ (sin afectar el resto) ---
p_subset_violin <- p_subset_violin +
  # intra: GTEX_AE (verde) y GTEX_AE_GEO (rojo)
  geom_segment(
    data = pair_intra,
    aes(x = x_br, xend = x_br, y = y_low, yend = y_high, color = col),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_segment(
    data = pair_intra,
    aes(x = x_br - tick, xend = x_br + tick, y = y_low, yend = y_low, color = col),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_segment(
    data = pair_intra,
    aes(x = x_br - tick, xend = x_br + tick, y = y_high, yend = y_high, color = col),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_text(
    data = pair_intra,
    aes(x = x_br, y = y_text, label = label, color = col),
    inherit.aes = FALSE,
    size = 4
  ) +
  # inter: AE_GEO vs GTEX_GEO (gris)
  geom_segment(
    data = pair_inter,
    aes(x = x_br, xend = x_br, y = y_low, yend = y_high, color = col),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_segment(
    data = pair_inter,
    aes(x = x_br - tick, xend = x_br + tick, y = y_low, yend = y_low, color = col),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_segment(
    data = pair_inter,
    aes(x = x_br - tick, xend = x_br + tick, y = y_high, yend = y_high, color = col),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_text(
    data = pair_inter,
    aes(x = x_br, y = y_text, label = label, color = col),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_color_identity()

print(p_subset_violin + coord_cartesian(ylim = c(0, 4000)))