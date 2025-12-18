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

adenosine_genes <- c(
  "ADORA1","ADORA2A","ADORA2B","ADORA3",
  "ENTPD1","NT5E",
  "ADA","ADK","PNP",
  "SLC29A1",
  "CD38","ENPP1"
)

padj_cut <- 0.05
lfc_cut  <- 1

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

pick_symbol_col <- function(df) {
  candidates <- c("Symbol", "Symbol_autho", "gene", "Gene")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) stop("No encuentro columna de símbolo (Symbol/Symbol_autho/Gene...).")
  hit[1]
}

prep_volcano_df <- function(df) {
  sym_col <- pick_symbol_col(df)
  
  df %>%
    dplyr::mutate(
      Gene = as.character(.data[[sym_col]]),
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", is.finite(log2FoldChange), !is.na(padj), padj > 0) %>%
    dplyr::mutate(
      neglog10_padj = -log10(padj),
      Type = dplyr::case_when(
        log2FoldChange >=  lfc_cut & padj < padj_cut ~ "Upregulated",
        log2FoldChange <= -lfc_cut & padj < padj_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      is_adenosine = Gene %in% adenosine_genes,
      is_adenosine_DE = is_adenosine & padj < padj_cut & abs(log2FoldChange) >= lfc_cut
    )
}

make_volcano <- function(res_df, title, subtitle, out_pdf) {
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = Type), size = 0.8, alpha = 0.8) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not Significant" = "gray70"
    )) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), color = "black", linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cut),     color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "log2FoldChange",
      y = expression(-log[10]("padj"))
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  ad_df <- res_df %>% dplyr::filter(is_adenosine_DE)
  
  if (nrow(ad_df) > 0) {
    p <- p +
      geom_point(
        data = ad_df,
        aes(x = log2FoldChange, y = neglog10_padj),
        inherit.aes = FALSE,
        color = "purple3",
        size = 1.4,
        alpha = 0.95
      ) +
      ggrepel::geom_text_repel(
        data = ad_df,
        aes(x = log2FoldChange, y = neglog10_padj, label = Gene),
        inherit.aes = FALSE,
        size = 3.2,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0
      )
  } else {
    message("Nota: No hay genes de adenosina que pasen |log2FC|>= ", lfc_cut,
            " y padj<", padj_cut, " para: ", subtitle)
  }
  
  ggsave(out_pdf, p, width = 8.5, height = 6.5, device = cairo_pdf)
  invisible(p)
}

# ===================== RUN ===================== #
DE_raw <- read_rds_safe(paths$raw)
DE_ae  <- read_rds_safe(paths$ae)

raw_df <- prep_volcano_df(DE_raw)
ae_df  <- prep_volcano_df(DE_ae)

make_volcano(
  res_df   = raw_df,
  title    = "Differential expression in ovarian cancer",
  subtitle = "Ovary tumors vs ovary control GTEx",
  out_pdf  = file.path(out_dir, "volcano_Ovary_control_GTEx_adenosine_highlight.pdf")
)

make_volcano(
  res_df   = ae_df,
  title    = "Differential expression in ovarian cancer",
  subtitle = "Ovary tumors vs reference control",
  out_pdf  = file.path(out_dir, "volcano_Reference_control_adenosine_highlight.pdf")
)

# Guardar tabla de genes de adenosina (aunque no pasen umbral), y marcar si pasan umbral
raw_ad <- raw_df %>%
  dplyr::filter(is_adenosine) %>%
  dplyr::select(Gene, log2FoldChange, padj, Type, is_adenosine_DE) %>%
  dplyr::arrange(dplyr::desc(is_adenosine_DE), Gene)

ae_ad <- ae_df %>%
  dplyr::filter(is_adenosine) %>%
  dplyr::select(Gene, log2FoldChange, padj, Type, is_adenosine_DE) %>%
  dplyr::arrange(dplyr::desc(is_adenosine_DE), Gene)

readr::write_tsv(raw_ad, file.path(out_dir, "adenosine_genes_status_Ovary_control_GTEx.tsv"))
readr::write_tsv(ae_ad,  file.path(out_dir, "adenosine_genes_status_Reference_control.tsv"))

message("\nDone. Outputs en: ", out_dir)
