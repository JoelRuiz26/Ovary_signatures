#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw = file.path(base_dir, "0_DGE_raw",         "DE_full_OVARY_DESeq2.rds"),
  ae  = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds")
)

out_dir <- file.path(base_dir, "2_2_DGE_volcano_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

adenosine_genes <- unique(c(
  # Receptores
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  
  # Transportadores (solo los mencionados en el documento)
  "SLC28A1",  # CNT
  "SLC29A1",  # ENT
  
  # Enzimas
  "NT5E",     # CD73
  "ENTPD1",   # CD39
  "DPP4",     # CD26
  "ADK",      # Adenosine kinase
  "ADA",      # Adenosine deaminase (ADA1 en el documento)
  "ADA2"      # Adenosine deaminase 2
))



padj_cut <- 0.05
lfc_cut  <- 1
max_neglog10 <- 301   # filtra -log10(p) > 301 antes del plot

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

make_volcano <- function(df, subtitle, out_pdf) {
  
  ad_df <- df %>%
    filter(is_adenosine_DE) %>%
    group_by(Gene) %>%
    slice_max(order_by = neglog10, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  p <- ggplot(df, aes(log2FoldChange, neglog10)) +
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
      aes(x = log2FoldChange, y = neglog10, label = Gene),  # <<< FIX
      inherit.aes = FALSE,
      size = 3.2,
      max.overlaps = Inf,
      fill = "white",
      color = "black",
      label.size = 0.25
    ) +
    labs(
      title = "Differential expression in ovarian cancer",
      subtitle = subtitle,
      x = "log2FoldChange",
      y = "-log10(padj)"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(out_pdf, p, width = 8, height = 6, device = cairo_pdf)
}

# ===================== RUN ===================== #
raw_df <- prep_df(read_rds_safe(paths$raw))
ae_df  <- prep_df(read_rds_safe(paths$ae))

make_volcano(
  raw_df,
  "Ovary tumors vs ovary control GTEx",
  file.path(out_dir, "volcano_Ovary_control_GTEx_adenosine.pdf")
)

make_volcano(
  ae_df,
  "Ovary tumors vs reference control",
  file.path(out_dir, "volcano_Reference_control_adenosine.pdf")
)

message("Done. Outputs en: ", out_dir)
