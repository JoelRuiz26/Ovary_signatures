#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(grid)
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

adenosine_genes <- c(
  "ADORA1","ADORA2A","ADORA2B","ADORA3",
  "ENTPD1","NT5E",
  "ADA","ADK","PNP",
  "SLC29A1",
  "CD38","ENPP1"
)

padj_cut <- 0.05
lfc_cut  <- 1
max_labels <- 30

# ===================== LOAD ===================== #
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe: ", p)
  readRDS(p)
}

prep_volcano_df <- function(df) {
  
  if (!("Symbol" %in% colnames(df))) {
    if ("Symbol_autho" %in% colnames(df)) {
      df <- df %>% mutate(Symbol = as.character(Symbol_autho))
    } else {
      stop("No encuentro Symbol ni Symbol_autho.")
    }
  }
  
  df %>%
    mutate(
      Gene           = as.character(Symbol),
      log2FoldChange = as.numeric(log2FoldChange),
      padj_used      = as.numeric(padj_safe),
      neglog10_padj  = as.numeric(mlog10_padj)
    ) %>%
    filter(
      !is.na(Gene),
      is.finite(log2FoldChange),
      is.finite(neglog10_padj),
      padj_used > 0
    ) %>%
    mutate(
      Type = case_when(
        log2FoldChange >=  lfc_cut & padj_used < padj_cut ~ "Upregulated",
        log2FoldChange <= -lfc_cut & padj_used < padj_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      is_adenosine = Gene %in% adenosine_genes,
      is_adenosine_DE = is_adenosine &
        padj_used < padj_cut &
        abs(log2FoldChange) >= lfc_cut
    )
}

# ===================== PLOT ===================== #
make_volcano <- function(df, subtitle, out_pdf) {
  
  ad_df <- df %>%
    filter(is_adenosine_DE) %>%
    arrange(desc(neglog10_padj)) %>%
    slice_head(n = max_labels)
  
  p <- ggplot(df, aes(log2FoldChange, neglog10_padj)) +
    geom_point(aes(color = Type), size = 0.8, alpha = 0.75) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not Significant" = "gray70"
    )) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cut),
               linetype = "dashed", linewidth = 0.4) +
    geom_point(
      data = ad_df,
      aes(log2FoldChange, neglog10_padj),
      inherit.aes = FALSE,
      shape = 21, fill = "purple3",
      color = "black", size = 2.8
    ) +
    ggrepel::geom_label_repel(
      data = ad_df,
      aes(label = Gene),
      size = 3.3,
      label.size = 0.25,
      box.padding = 0.35,
      point.padding = 0.2,
      max.overlaps = Inf
    ) +
    coord_cartesian(ylim = c(0, 300)) +   # <<< ÚNICA CORRECCIÓN
    labs(
      title = "Differential expression in ovarian cancer",
      subtitle = subtitle,
      x = "log2FoldChange",
      y = expression(-log[10]("padj_safe"))
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(out_pdf, p, width = 8.5, height = 6.5, device = cairo_pdf)
}

# ===================== RUN ===================== #
raw_df <- prep_volcano_df(read_rds_safe(paths$raw))
ae_df  <- prep_volcano_df(read_rds_safe(paths$ae))

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

# ===================== SAVE TABLES ===================== #
raw_df %>%
  filter(is_adenosine) %>%
  select(Gene, log2FoldChange, padj_used, neglog10_padj, Type, is_adenosine_DE) %>%
  arrange(desc(is_adenosine_DE), Gene) %>%
  write_tsv(file.path(out_dir, "adenosine_genes_status_Ovary_control_GTEx.tsv"))

ae_df %>%
  filter(is_adenosine) %>%
  select(Gene, log2FoldChange, padj_used, neglog10_padj, Type, is_adenosine_DE) %>%
  arrange(desc(is_adenosine_DE), Gene) %>%
  write_tsv(file.path(out_dir, "adenosine_genes_status_Reference_control.tsv"))

message("Done.")
