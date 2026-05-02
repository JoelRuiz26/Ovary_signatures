#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(ggplot2)
  library(tibble)
})

IN_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_Counts_metadata/"
OUT_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
# =========================
# 1) Load normalized matrices
# =========================
norm_all_ovary <- readRDS(file.path(IN_DIR, "3_3_2_norm_counts_tumor_GTEX.rds"))
norm_all_auto  <- readRDS(file.path(IN_DIR, "3_3_2_norm_count_tumor_AE.rds"))

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

print(summary(fano_AE$fano_factor))
print(summary(fano_GTEX$fano_factor))

# =========================
# 5) Identify gene subsets (by presence flags ae/gtex/geo)
# =========================
all_genes <- vroom::vroom("~/Ovary_signatures/3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv")

# normaliza nombres y arregla duplicados (ae/gtex/geo salen 2 veces: dir + 0/1)
names(all_genes) <- tolower(names(all_genes))
names(all_genes) <- make.unique(names(all_genes))

# arregla duplicados: primero renombra columnas de dirección
all_genes <- all_genes %>%
  dplyr::rename(
    ae_dir   = ae,
    gtex_dir = gtex,
    geo_dir  = geo
  ) %>%
  dplyr::rename(
    ae   = ae.1,
    gtex = gtex.1,
    geo  = geo.1
  ) %>%
  dplyr::mutate(
    ae   = as.integer(ae),
    gtex = as.integer(gtex),
    geo  = as.integer(geo)
  )

# etiqueta de subconjunto por presencia/ausencia
map_src <- all_genes %>%
  dplyr::transmute(
    ensembl = ensembl,
    gene = gene,
    dataset = dplyr::case_when(
      ae == 1 & gtex == 0 & geo == 0 ~ "AE",
      ae == 0 & gtex == 1 & geo == 0 ~ "GTEX",
      ae == 0 & gtex == 0 & geo == 1 ~ "GEO",
      ae == 1 & gtex == 1 & geo == 0 ~ "GTEX_AE",
      ae == 1 & gtex == 0 & geo == 1 ~ "AE_GEO",
      ae == 0 & gtex == 1 & geo == 1 ~ "GTEX_GEO",
      ae == 1 & gtex == 1 & geo == 1 ~ "GTEX_AE_GEO",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::distinct(ensembl, dataset)

# pega el subconjunto a ambos fano
fano_GTEX <- fano_GTEX %>% dplyr::left_join(map_src, by = c("gene" = "ensembl"))
fano_AE   <- fano_AE   %>% dplyr::left_join(map_src, by = c("gene" = "ensembl"))
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
  "GTEX_AE_GEO" = "#F9B8B9"
)


# =========================
# 7) PLOTS (LINEAR axis)
#    IMPORTANT CHANGE:
#    When list == "AE", include all genes with ae==1 (i.e., any dataset label containing AE)
#    When list == "GTEx", include all genes with gtex==1 (any dataset label containing GTEX)
#    IMPORTANT for deltas:
#    - Replace the "AE vs GTEX" separated violin to use TOTALS (AE_total vs GTEX_total)
#      so the Δ updates according to your new definition.
# =========================

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

# =========================
# 7.3) VIOLIN (tu violín original + Δ)
# =========================
df_violin <- df_subset_zoom %>%
  mutate(
    dataset = factor(
      dataset,
      levels = c("AE","GTEX","AE_GEO","GTEX_GEO","GTEX_AE","GTEX_AE_GEO")
    ),
    list = factor(list, levels = c("GTEx", "AE"))
  )

DODGE_W <- 1

# --- plot base (NO cambies estética, solo el dodge) ---
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
    position  = position_dodge(width = DODGE_W)
  ) +
  geom_boxplot(
    width         = 0.12,
    outlier.shape = NA,
    linewidth     = 0.25,
    fill          = NA,
    position      = position_dodge(width = DODGE_W)
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
# Brackets + Δ (como lo tenías)
# =========================
x_base <- setNames(seq_along(levels(df_violin$dataset)), levels(df_violin$dataset))

med_intra <- df_violin %>%
  filter(dataset %in% c("GTEX_AE", "GTEX_AE_GEO")) %>%
  group_by(dataset, list) %>%
  summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
  mutate(x0 = x_base[as.character(dataset)])

offset <- DODGE_W / 4

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
    col = "gray30",
    x_br = x0,
    tick = 0.10,
    y_text = y_high + 400,
    label = paste0("\u0394=", round(delta, 1))
  )

med_inter <- df_violin %>%
  filter(dataset %in% c("AE_GEO", "GTEX_GEO")) %>%
  group_by(dataset) %>%
  summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
  mutate(x = x_base[as.character(dataset)])

pair_inter <- tibble(
  x_br  = mean(med_inter$x),
  y_low = min(med_inter$med),
  y_high= max(med_inter$med),
  delta = max(med_inter$med) - min(med_inter$med),
  col   = "gray30",
  tick  = 0.12,
  y_text = max(med_inter$med) + 300,
  label = paste0("\u0394=", round(delta, 1))
)

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

print(p_subset_violin + coord_cartesian(ylim = c(0, 12000)))

# =========================================================
# 7.4) NUEVOS VIOLINES SEPARADOS (MISMA estética)
#    MINIMAL CHANGE:
#    - Panel (1) is now TOTALS: AE_total vs GTEX_total (so Δ updates)
#    - Panels (2)-(4) stay identical
# =========================================================

make_base_violin <- function(df_in, title_txt, ylim_max = 12000) {
  ggplot(
    df_in,
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
      position  = position_dodge(width = DODGE_W)
    ) +
    geom_boxplot(
      width         = 0.12,
      outlier.shape = NA,
      linewidth     = 0.25,
      fill          = NA,
      position      = position_dodge(width = DODGE_W)
    ) +
    scale_fill_manual(values = palette_dataset, na.value = "grey80") +
    labs(
      title = title_txt,
      x = "Gene subset",
      y = "Fano factor"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0, ylim_max))
}

add_intra_delta_single <- function(p, df_in) {
  x_base_local <- setNames(seq_along(levels(df_in$dataset)), levels(df_in$dataset))
  
  med_intra_local <- df_in %>%
    group_by(dataset, list) %>%
    summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
    mutate(x0 = x_base_local[as.character(dataset)])
  
  offset_local <- DODGE_W / 4
  
  pair_intra_local <- med_intra_local %>%
    mutate(x = ifelse(list == "GTEx", x0 - offset_local, x0 + offset_local)) %>%
    group_by(dataset, x0) %>%
    summarise(
      y_low  = min(med),
      y_high = max(med),
      delta  = y_high - y_low,
      .groups = "drop"
    ) %>%
    mutate(
      col = "gray30",
      x_br = x0,
      tick = 0.05,
      y_text = y_high + 400,
      label = paste0("\u0394=", round(delta, 1))
    )
  
  p +
    geom_segment(
      data = pair_intra_local,
      aes(x = x_br, xend = x_br, y = y_low, yend = y_high, color = col),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_segment(
      data = pair_intra_local,
      aes(x = x_br - tick, xend = x_br + tick, y = y_low, yend = y_low, color = col),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_segment(
      data = pair_intra_local,
      aes(x = x_br - tick, xend = x_br + tick, y = y_high, yend = y_high, color = col),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_text(
      data = pair_intra_local,
      aes(x = x_br, y = y_text, label = label, color = col),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_color_identity()
}

add_inter_delta_pair <- function(p, df_in) {
  x_base_local <- setNames(seq_along(levels(df_in$dataset)), levels(df_in$dataset))
  
  med_inter_local <- df_in %>%
    group_by(dataset) %>%
    summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
    mutate(x = x_base_local[as.character(dataset)])
  
  pair_inter_local <- tibble(
    x_br   = mean(med_inter_local$x),
    y_low  = min(med_inter_local$med),
    y_high = max(med_inter_local$med),
    delta  = max(med_inter_local$med) - min(med_inter_local$med),
    col    = "gray30",
    tick   = 0.12,
    y_text = max(med_inter_local$med) + 300,
    label  = paste0("\u0394=", round(delta, 1))
  )
  
  p +
    geom_segment(
      data = pair_inter_local,
      aes(x = x_br, xend = x_br, y = y_low, yend = y_high, color = col),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_segment(
      data = pair_inter_local,
      aes(x = x_br - tick, xend = x_br + tick, y = y_low, yend = y_low, color = col),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_segment(
      data = pair_inter_local,
      aes(x = x_br - tick, xend = x_br + tick, y = y_high, yend = y_high, color = col),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    geom_text(
      data = pair_inter_local,
      aes(x = x_br, y = y_text, label = label, color = col),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_color_identity()
}

# ---- (1) AE_total vs GTEX_total (REPLACES old AE vs GTEX exclusive panel) ----
df_v1 <- df_subset_zoom %>%   # <-- antes era df_plot_subset
  mutate(
    dataset = ifelse(list == "AE", "AE_total", "GTEX_total"),
    dataset = factor(dataset, levels = c("AE_total", "GTEX_total"))
  ) %>%
  droplevels()

# same aesthetics but with a small local palette for totals
palette_totals <- c("AE_total" = "#F4D96B", "GTEX_total" = "#7A8DB8")

p_violin_AE_GTEX <- ggplot(
  df_v1,
  aes(
    x = dataset,
    y = fano_factor,
    fill = dataset
  )
) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = NA) +
  scale_fill_manual(values = palette_totals, na.value = "grey80") +
  labs(
    title = "Fano factor (AE_total vs GTEX_total)",
    x = "Gene subset",
    y = "Fano factor"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 12000))

# Δ between totals
med_tot <- df_v1 %>%
  group_by(dataset) %>%
  summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop") %>%
  mutate(x = as.numeric(dataset))

pair_tot <- tibble(
  x_br   = mean(med_tot$x),
  y_low  = min(med_tot$med),
  y_high = max(med_tot$med),
  delta  = max(med_tot$med) - min(med_tot$med),
  col    = "gray30",
  tick   = 0.05,
  y_text = max(med_tot$med) + 300,
  label  = paste0("\u0394=", round(delta, 1))
)

p_violin_AE_GTEX <- p_violin_AE_GTEX +
  geom_segment(
    data = pair_tot,
    aes(x = x_br, xend = x_br, y = y_low, yend = y_high),
    inherit.aes = FALSE,
    linewidth = 0.9,
    color = "gray30"
  ) +
  geom_segment(
    data = pair_tot,
    aes(x = x_br - tick, xend = x_br + tick, y = y_low, yend = y_low),
    inherit.aes = FALSE,
    linewidth = 0.9,
    color = "gray30"
  ) +
  geom_segment(
    data = pair_tot,
    aes(x = x_br - tick, xend = x_br + tick, y = y_high, yend = y_high),
    inherit.aes = FALSE,
    linewidth = 0.9,
    color = "gray30"
  ) +
  geom_text(
    data = pair_tot,
    aes(x = x_br, y = y_text, label = label),
    inherit.aes = FALSE,
    size = 4,
    color = "gray30"
  )

print(p_violin_AE_GTEX)

# ---- (2) GTEX_AE (verde) ----
df_v2 <- df_violin %>%
  filter(dataset == "GTEX_AE") %>%
  droplevels()

p_violin_GTEX_AE <- make_base_violin(
  df_in = df_v2,
  title_txt = "Fano factor (GTEX_AE)",
  ylim_max = 12000
)

if (all(c("GTEx", "AE") %in% levels(df_v2$list))) {
  p_violin_GTEX_AE <- add_intra_delta_single(p_violin_GTEX_AE, df_v2)
}
print(p_violin_GTEX_AE)

# ---- (3) GTEX_AE_GEO (rojo) ----
df_v3 <- df_violin %>%
  filter(dataset == "GTEX_AE_GEO") %>%
  droplevels()

p_violin_GTEX_AE_GEO <- make_base_violin(
  df_in = df_v3,
  title_txt = "Fano factor (GTEX_AE_GEO)",
  ylim_max = 12000
)

if (all(c("GTEx", "AE") %in% levels(df_v3$list))) {
  p_violin_GTEX_AE_GEO <- add_intra_delta_single(p_violin_GTEX_AE_GEO, df_v3)
}
print(p_violin_GTEX_AE_GEO)

# ---- (4) AE_GEO vs GTEX_GEO ----
df_v4 <- df_violin %>%
  filter(dataset %in% c("AE_GEO", "GTEX_GEO")) %>%
  mutate(dataset = factor(dataset, levels = c("AE_GEO", "GTEX_GEO"))) %>%
  droplevels()

p_violin_AE_GEO_vs_GTEX_GEO <- make_base_violin(
  df_in = df_v4,
  title_txt = "Fano factor (AE_GEO vs GTEX_GEO)",
  ylim_max = 12000
)

if (all(c("AE_GEO", "GTEX_GEO") %in% levels(df_v4$dataset))) {
  p_violin_AE_GEO_vs_GTEX_GEO <- add_inter_delta_pair(p_violin_AE_GEO_vs_GTEX_GEO, df_v4)
}
print(p_violin_AE_GEO_vs_GTEX_GEO)

# =========================
# 8) SAVE PLOTS
# =========================
save_plot_both <- function(p, file_stub, width = 9, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, paste0(file_stub, ".png")),
    plot     = p,
    width    = width,
    height   = height,
    units    = "in",
    dpi      = dpi
  )
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, paste0(file_stub, ".pdf")),
    plot     = p,
    width    = width,
    height   = height,
    units    = "in",
    device   = cairo_pdf
  )
}

# ---- Guardar LINES ----
p_subset_0_4000 <- p_subset + coord_cartesian(xlim = c(0, 4000))
save_plot_both(
  p         = p_subset_0_4000,
  file_stub = "3_3_4_fano_subset_lines",
  width     = 10,
  height    = 6
)

# ---- Guardar VIOLIN original con Δ ----
p_violin_0_14000 <- p_subset_violin + coord_cartesian(ylim = c(0, 12000))
save_plot_both(
  p         = p_violin_0_14000,
  file_stub = "3_3_4_fano_subset_violin_delta",
  width     = 11,
  height    = 6
)

# ---- Guardar NUEVOS violines separados ----
save_plot_both(
  p         = p_violin_AE_GTEX,
  file_stub = "3_3_4_fano_violin_AE_total_vs_GTEX_total",
  width     = 8,
  height    = 6
)

save_plot_both(
  p         = p_violin_GTEX_AE,
  file_stub = "3_3_4_fano_violin_GTEX_AE_green",
  width     = 6,
  height    = 6
)

save_plot_both(
  p         = p_violin_GTEX_AE_GEO,
  file_stub = "3_3_4_fano_violin_GTEX_AE_GEO_red",
  width     = 6,
  height    = 6
)

save_plot_both(
  p         = p_violin_AE_GEO_vs_GTEX_GEO,
  file_stub = "3_3_4_fano_violin_AE_GEO_vs_GTEX_GEO",
  width     = 8,
  height    = 6
)




