#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(patchwork)   # <-- grid + shared legend
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  # GEO: lista con 6 datasets
  geo = file.path(base_dir, "2_DEG_GEO",  "DEGs_alldsets_GEO.rds")
)

out_dir <- file.path(base_dir, "3_Consensus_DGE_analysis/3_2_DGE_volcano_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

adenosine_genes <- unique(c(
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  "SLC28A1", "SLC29A1",
  "NT5E", "ENTPD1", "DPP4", "ADK", "ADA", "ADA2"
))

padj_cut <- 0.05
lfc_cut  <- 1
max_neglog10 <- 301   # filtra -log10(p) > 301 antes del plot

p_floor <- 1e-300      # evitar log(0) en -log10

# ===================== HELPERS ===================== #
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe: ", p)
  readRDS(p)
}

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
      neglog10 = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10) %>%
    mutate(
      Type = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Upregulated",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      is_adenosine_DE =
        Gene %in% adenosine_genes &
        padj_used < padj_cut &
        abs(log2FoldChange) >= lfc_cut
    )
}

# ===== GEO (limma) =====
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
      neglog10 = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10)
  
  # múltiples probes por gen → deja 1 fila por Gene (la de mayor |logFC|)
  out <- out %>%
    group_by(Gene) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Type = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Upregulated",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      is_adenosine_DE =
        Gene %in% adenosine_genes &
        padj_used < padj_cut &
        abs(log2FoldChange) >= lfc_cut
    )
  
  out
}

# ---- Make ONE volcano plot (returns ggplot object) ----
make_volcano_plot <- function(df, subtitle) {
  
  ad_df <- df %>%
    filter(is_adenosine_DE) %>%
    group_by(Gene) %>%
    slice_max(order_by = neglog10, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ggplot(df, aes(log2FoldChange, neglog10)) +
    geom_point(aes(color = Type), size = 0.8, alpha = 0.75) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not Significant" = "gray70"
    )) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cut),     linetype = "dashed", linewidth = 0.4) +
    geom_point(
      data = ad_df,
      aes(x = log2FoldChange, y = neglog10),
      inherit.aes = FALSE,
      shape = 21,
      fill  = "purple3",
      color = "black",
      stroke = 0.35,
      size = 2.8,
      alpha = 0.98
    ) +
    ggrepel::geom_label_repel(
      data = ad_df,
      aes(x = log2FoldChange, y = neglog10, label = Gene),
      inherit.aes = FALSE,
      size = 3.2,
      max.overlaps = Inf,
      fill = "white",
      color = "black",
      label.size = 0.25
    ) +
    labs(
      title = subtitle,
      x = "log2FoldChange",
      y = "-log10(padj)",
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    )
}

# ---- Save individual volcano (PDF+PNG) ----
save_volcano <- function(p, out_pdf, w = 8, h = 6) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_pdf, p, width = w, height = h, device = cairo_pdf)
  out_png <- sub("\\.pdf$", ".png", out_pdf)
  ggsave(out_png, p, width = w, height = h, dpi = 600, device = "png", bg = "white")
}

# ===================== RUN ===================== #
plots <- list()

# ---- GTEx ----
raw_df <- prep_df(read_rds_safe(paths$raw))
p_gtex <- make_volcano_plot(raw_df, "GTEx control")
save_volcano(p_gtex, file.path(out_dir, "3_2_1_volcano_Ovary_control_GTEx.pdf"))
plots[["GTEx"]] <- p_gtex

# ---- AE ----
ae_df <- prep_df(read_rds_safe(paths$ae))
p_ae <- make_volcano_plot(ae_df, "Reference control (AE)")
save_volcano(p_ae, file.path(out_dir, "3_2_1_volcano_Reference_control_AE.pdf"))
plots[["AE"]] <- p_ae

# ---- GEO (6 volcanos) ----
GEO_full <- NULL
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
    message("GEO volcano: ", nm)
    
    df_geo <- prep_df_geo(GEO_full[[nm]])
    nm_safe <- gsub("[^A-Za-z0-9_\\-]", "_", nm)
    
    p_geo <- make_volcano_plot(df_geo, nm_safe)
    save_volcano(p_geo, file.path(out_dir, paste0("3_2_1_volcano_GEO_", nm_safe, ".pdf")))
    
    plots[[paste0("GEO_", nm_safe)]] <- p_geo
  }
} else {
  message("Nota: no existe GEO RDS en: ", paths$geo, " (se omiten volcanos GEO).")
}

# ---- GRID (todos + una sola leyenda) ----
# Nota: si por alguna razón no son 8 (p.ej., GEO no cargó), igual arma el grid con lo que haya.
grid <- wrap_plots(plots, ncol = 4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

grid_pdf <- file.path(out_dir, "3_2_1_volcano_GRID_GTEx_AE_GEO.pdf")
grid_png <- sub("\\.pdf$", ".png", grid_pdf)

ggsave(grid_pdf, grid, width = 16, height = 8.5, units = "in", device = cairo_pdf)
ggsave(grid_png, grid, width = 16, height = 8.5, units = "in", dpi = 600, device = "png", bg = "white")

message("Done. Outputs en: ", out_dir)
