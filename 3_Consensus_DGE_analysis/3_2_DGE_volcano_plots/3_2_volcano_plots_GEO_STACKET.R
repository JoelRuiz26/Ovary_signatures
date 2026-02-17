#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(patchwork)   # grid + shared legend
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  geo = file.path(base_dir, "2_DEG_GEO",  "DEGs_alldsets_GEO.rds")  # lista con 6 datasets
)

out_dir <- file.path(base_dir, "3_Consensus_DGE_analysis/3_2_DGE_volcano_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cut <- 0.01
lfc_cut  <- 1
TOP_N    <- 5

max_neglog10 <- 301   # filtra -log10(p) > 301 antes del plot
p_floor <- 1e-300     # evitar log(0) en -log10

# ===================== HELPERS ===================== #
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe: ", p)
  readRDS(p)
}

# ---- DESeq2-like prep ----
prep_df <- function(df) {
  
  if (!("Symbol" %in% colnames(df))) {
    if ("Symbol_autho" %in% colnames(df)) {
      df <- df %>% mutate(Symbol = as.character(Symbol_autho))
    } else {
      stop("No encuentro Symbol ni Symbol_autho.")
    }
  }
  
  padj_col <- if ("padj_safe" %in% colnames(df)) "padj_safe" else if ("padj" %in% colnames(df)) "padj" else NA
  if (is.na(padj_col)) stop("No encuentro padj_safe ni padj.")
  
  df %>%
    transmute(
      Gene = as.character(Symbol),
      log2FoldChange = as.numeric(log2FoldChange),
      padj_used = as.numeric(.data[[padj_col]])
    ) %>%
    filter(
      !is.na(Gene), Gene != "",
      is.finite(log2FoldChange),
      is.finite(padj_used),
      padj_used > 0
    ) %>%
    mutate(
      padj_used = pmax(padj_used, p_floor),
      neglog10  = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10) %>%
    mutate(
      Type = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Overexpressed",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Underexpressed",
        TRUE ~ "Not significant"
      ),
      direction = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Up",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ NA_character_
      ),
      score = abs(log2FoldChange) * neglog10
    )
}

# ---- GEO (limma-like prep) ----
prep_df_geo <- function(df) {
  
  gene_col <- dplyr::case_when(
    "Symbol" %in% colnames(df)       ~ "Symbol",
    "Symbol_autho" %in% colnames(df) ~ "Symbol_autho",
    "Gene" %in% colnames(df)         ~ "Gene",
    "gene" %in% colnames(df)         ~ "gene",
    "Gene.symbol" %in% colnames(df)  ~ "Gene.symbol",
    "genesymbol" %in% colnames(df)   ~ "genesymbol",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col)) stop("GEO: no encuentro columna de símbolo (Symbol/Symbol_autho/Gene/Gene.symbol/...).")
  
  lfc_col <- if ("logFC" %in% colnames(df)) "logFC" else if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else NA
  if (is.na(lfc_col)) stop("GEO: no encuentro columna de logFC (logFC/log2FoldChange).")
  
  padj_col <- dplyr::case_when(
    "adj.P.Val" %in% colnames(df) ~ "adj.P.Val",
    "FDR" %in% colnames(df)       ~ "FDR",
    "padj" %in% colnames(df)      ~ "padj",
    "padj_safe" %in% colnames(df) ~ "padj_safe",
    TRUE ~ NA_character_
  )
  if (is.na(padj_col)) stop("GEO: no encuentro columna de p ajustada (adj.P.Val/FDR/padj/...).")
  
  out <- df %>%
    transmute(
      Gene = as.character(.data[[gene_col]]),
      log2FoldChange = as.numeric(.data[[lfc_col]]),
      padj_used = as.numeric(.data[[padj_col]])
    ) %>%
    filter(
      !is.na(Gene), Gene != "",
      is.finite(log2FoldChange),
      is.finite(padj_used),
      padj_used > 0
    ) %>%
    mutate(
      padj_used = pmax(padj_used, p_floor),
      neglog10  = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10)
  
  # múltiples probes por gen → deja 1 fila por Gene (la de mayor |logFC|)
  out <- out %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Type = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Overexpressed",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Underexpressed",
        TRUE ~ "Not significant"
      ),
      direction = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Up",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ NA_character_
      ),
      score = abs(log2FoldChange) * neglog10
    )
  
  out
}

# ---- pick Top genes to label: Top N Up + Top N Down by score ----
pick_top_labels <- function(df, top_n = 5) {
  
  cand <- df %>%
    filter(!is.na(direction)) %>%
    filter(!grepl("-", Gene) & !grepl("\\.", Gene)) %>%   # evita guion y punto
    group_by(Gene) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  top_up <- cand %>%
    filter(direction == "Up") %>%
    slice_max(order_by = score, n = top_n, with_ties = FALSE)
  
  top_down <- cand %>%
    filter(direction == "Down") %>%
    slice_max(order_by = score, n = top_n, with_ties = FALSE)
  
  bind_rows(top_up, top_down)
}

# ===================== RUN ===================== #
# ---- build ONE long df: all datasets stacked ----
all_df <- tibble()

# GTEx
raw_df <- prep_df(read_rds_safe(paths$raw)) %>% mutate(dataset = "GTEx control")
all_df <- bind_rows(all_df, raw_df)

# AE
ae_df <- prep_df(read_rds_safe(paths$ae)) %>% mutate(dataset = "Reference control (AE)")
all_df <- bind_rows(all_df, ae_df)

# GEO
if (file.exists(paths$geo)) {
  GEO_full <- read_rds_safe(paths$geo)
  
  if (is.data.frame(GEO_full)) {
    GEO_full <- list(GEO = GEO_full)
  }
  
  if (!is.list(GEO_full) || length(GEO_full) == 0) {
    stop("El objeto GEO no es una lista no-vacía. Revisa: ", paths$geo)
  }
  
  geo_names <- names(GEO_full)
  if (is.null(geo_names) || any(geo_names == "")) {
    geo_names <- paste0("GEO_", seq_along(GEO_full))
  }
  
  for (nm in geo_names) {
    message("GEO add: ", nm)
    nm_safe <- gsub("[^A-Za-z0-9_\\-]", "_", nm)
    df_geo  <- prep_df_geo(GEO_full[[nm]]) %>% mutate(dataset = nm_safe)
    all_df  <- bind_rows(all_df, df_geo)
  }
} else {
  message("Nota: no existe GEO RDS en: ", paths$geo, " (se omite GEO).")
}

# ---- consistent facet order (GTEx, AE, luego GEO) ----
dataset_levels <- c(
  setdiff(unique(all_df$dataset), c("GTEx control", "Reference control (AE)")),
  "GTEx control",
  "Reference control (AE)"
)

all_df <- all_df %>% mutate(dataset = factor(dataset, levels = dataset_levels))

# ---- global X axis for all panels ----
x_max <- max(abs(all_df$log2FoldChange), na.rm = TRUE)
x_lim <- c(-x_max, x_max)

# ---- label df: top genes per dataset ----
label_df <- all_df %>%
  group_by(dataset) %>%
  group_modify(~ pick_top_labels(.x, top_n = TOP_N)) %>%
  ungroup()

# ===================== PLOT: VOLCANOS EN FILA (dataset ~ .) ===================== #
p <- ggplot(all_df, aes(log2FoldChange, neglog10)) +
  geom_point(aes(color = Type), size = 0.8, alpha = 0.75) +
  scale_color_manual(
    values = c(
      "Overexpressed"   = "#D7301F",
      "Underexpressed"  = "#2C7FB8",
      "Not significant" = "grey75"
    ),
    breaks = c("Overexpressed", "Not significant", "Underexpressed")
  ) +
  guides(color = guide_legend(
    override.aes = list(shape = 16, size = 4, alpha = 1)
  )) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut),     linetype = "dashed", linewidth = 0.4) +
  geom_point(
    data = label_df,
    aes(x = log2FoldChange, y = neglog10, color = Type),
    inherit.aes = FALSE,
    size = 1.6,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  ggrepel::geom_label_repel(
    data = label_df,
    aes(x = log2FoldChange, y = neglog10, label = Gene),
    inherit.aes = FALSE,
    size = 3.2,
    max.overlaps = Inf,
    fill = "white",
    color = "black",
    label.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 1
  ) +
  facet_grid(dataset ~ ., scales = "free_y") +   # <-- CAMBIO CLAVE
  coord_cartesian(xlim = x_lim, clip = "off") +
  labs(
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "Expression"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text.y = element_text(face = "bold", size = 12),
    plot.margin = margin(10, 35, 10, 35)
  )

# ---- SAVE (una sola figura apilada) ----
out_pdf <- file.path(out_dir, "3_2_1_volcano_STACK_byDataset.pdf")
out_png <- sub("\\.pdf$", ".png", out_pdf)

# altura dinámica (≈ 2.3 in por panel)
n_panels <- length(levels(all_df$dataset))
ggsave(out_pdf, p, width = 8.5, height = max(6, 2.3 * n_panels), units = "in", device = cairo_pdf)
ggsave(out_png, p, width = 8.5, height = max(6, 2.3 * n_panels), units = "in", dpi = 600, device = "png", bg = "white")

message("Done. Outputs en: ", out_dir)
