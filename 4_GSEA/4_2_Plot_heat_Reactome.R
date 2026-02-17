#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ============================================================
# CONFIG
# ============================================================
base_dir <- path.expand("~/Ovary_signatures/4_GSEA/4_0_REACTOME")
in_dir   <- base_dir
out_dir  <- file.path(base_dir, "4_1_Heatmap_Top5UpDown")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pattern_in <- "_REACTOME_FULL\\.tsv$"

TOP_N_PER_SIDE <- 5
PADJ_CUT       <- 0.05

out_tbl_top  <- file.path(out_dir, "4_1_Top5UpDown_byDataset.tsv")
out_rds_plot <- file.path(out_dir, "4_2_Heatmap_REACTOME_Top5UpDown_plot.rds")
out_png      <- file.path(out_dir, "4_2_Heatmap_REACTOME_Top5UpDown.png")
out_pdf      <- file.path(out_dir, "4_2_Heatmap_REACTOME_Top5UpDown.pdf")

# ============================================================
# 1) READ + STANDARDIZE
# ============================================================
files <- list.files(in_dir, pattern = pattern_in, full.names = TRUE)
if (length(files) == 0) {
  stop("No encontré archivos que coincidan con: ", pattern_in, "\nRevisa in_dir = ", in_dir)
}

read_one <- function(path) {
  dataset <- basename(path) %>%
    str_replace(pattern_in, "") %>%
    str_replace("_REACTOME_FULL$", "")
  
  df <- read_tsv(path, col_types = cols(), progress = FALSE)
  
  required <- c("Description", "NES", "p.adjust")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(
      "Archivo: ", basename(path),
      "\nNo tiene columnas requeridas: ", paste(missing, collapse = ", "),
      "\nColumnas disponibles: ", paste(colnames(df), collapse = ", ")
    )
  }
  
  df %>%
    transmute(
      dataset     = dataset,
      Description = as.character(Description),
      NES         = as.numeric(NES),
      p.adjust    = as.numeric(p.adjust)
    ) %>%
    filter(!is.na(Description), Description != "",
           !is.na(NES), is.finite(NES),
           !is.na(p.adjust), is.finite(p.adjust)) %>%
    mutate(
      term_simple = Description %>%
        str_remove("^REACTOME_") %>%
        str_replace_all("_", " "),
      direction = ifelse(NES > 0, "Up", "Down")
    )
}

gsea_full <- map_dfr(files, read_one)

# ============================================================
# 2) FILTER SIGNIFICANT + TOP5 UP/DOWN PER DATASET
# ============================================================
x <- gsea_full %>%
  filter(p.adjust < PADJ_CUT)

if (nrow(x) == 0) {
  stop("Después de filtrar p.adjust < ", PADJ_CUT, " no quedó ningún término.")
}

top5_up_down_by_dataset <- x %>%
  group_by(dataset, direction) %>%
  arrange(desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = TOP_N_PER_SIDE) %>%
  ungroup()

write_tsv(top5_up_down_by_dataset, out_tbl_top)

# ============================================================
# 3) PREP FOR PLOT (orden de términos + estrellas + alpha)
# ============================================================
# Orden de datasets (AE, GTEx, luego GEOs)
geo_ds <- sort(unique(top5_up_down_by_dataset$dataset)[grepl("^GEO_", unique(top5_up_down_by_dataset$dataset))])

dataset_levels <- c("AE", "GTEx", geo_ds)
dataset_levels <- dataset_levels[dataset_levels %in% unique(top5_up_down_by_dataset$dataset)]

plot_df <- top5_up_down_by_dataset %>%
  group_by(term_simple) %>%
  mutate(absNES_max = max(abs(NES), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    dataset_chr = as.character(dataset),
    dataset     = factor(dataset_chr, levels = dataset_levels),
    sig         = -log10(p.adjust),
    
    # estrellas (muy legible)
    star = case_when(
      p.adjust <= 1e-3 ~ "***",
      p.adjust <= 1e-2 ~ "**",
      p.adjust <= 5e-2 ~ "*",
      TRUE             ~ ""
    )
  )

# alpha continuo: rescale global
sig_rng <- range(plot_df$sig, finite = TRUE)
if (!all(is.finite(sig_rng)) || diff(sig_rng) == 0) {
  plot_df <- plot_df %>% mutate(alpha_sig = 1.0)
} else {
  plot_df <- plot_df %>% mutate(alpha_sig = scales::rescale(sig, to = c(0.25, 1.0), from = sig_rng))
}

# Orden del eje Y: por absNES_max (estable)
plot_df <- plot_df %>%
  mutate(term_simple = reorder(term_simple, absNES_max))

# ✅ Labels para facet (evita solape): GEO\nGSE#####
plot_df <- plot_df %>%
  mutate(
    dataset_label = case_when(
      dataset_chr %in% c("AE", "GTEx") ~ dataset_chr,
      TRUE ~ str_replace(dataset_chr, "^(GEO)_", "\\1\n")
    )
  )

# niveles del label en el mismo orden que dataset_levels
dataset_label_levels <- c(
  intersect(c("AE", "GTEx"), dataset_levels),
  str_replace(geo_ds, "^(GEO)_", "\\1\n")
)
dataset_label_levels <- dataset_label_levels[dataset_label_levels %in% unique(plot_df$dataset_label)]

plot_df <- plot_df %>%
  mutate(dataset_label = factor(dataset_label, levels = dataset_label_levels))

# ============================================================
# 4) PLOT (estilo igual al primero)
# ============================================================
p <- ggplot(plot_df, aes(dataset_label, term_simple)) +
  geom_tile(
    aes(fill = NES, alpha = alpha_sig),
    color = "grey92", linewidth = 0.4
  ) +
  geom_text(
    aes(label = star),
    fontface = "bold", size = 4, color = "grey15"
  ) +
  
  # 🎨 PALETA ELEGANTE (azul-blanco-rojo)
  scale_fill_gradient2(
    name = "NES",
    low  = "#3B4CC0",
    mid  = "#F7F7F7",
    high = "#B40426",
    midpoint = 0
  ) +
  
  scale_alpha(guide = "none") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    
    axis.text.y = element_text(size = 10, face = "bold"),
    
    # quitamos eje X (usamos strip)
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold"),
    
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(
      size = 12,
      face = "bold",
      lineheight = 0.95,
      margin = margin(t = 6, b = 6)
    ),
    
    plot.margin = margin(t = 10, r = 20, b = 18, l = 10)
  ) +
  facet_grid(. ~ dataset_label, switch = "x")

# ============================================================
# 5) SAVE
# ============================================================
saveRDS(p, out_rds_plot)

n_ds <- length(unique(plot_df$dataset_label))

ggsave(
  out_png,
  p,
  width  = max(14, 2.2 * n_ds),
  height = 9,
  units  = "in",
  dpi    = 600,
  limitsize = FALSE
)

ggsave(
  out_pdf,
  p,
  width  = max(14, 2.2 * n_ds),
  height = 9,
  units  = "in",
  limitsize = FALSE,
  device = cairo_pdf
)

message("\nDONE ✅")
message("Top table: ", out_tbl_top)
message("Plot RDS : ", out_rds_plot)
message("PNG      : ", out_png)
message("PDF      : ", out_pdf)
