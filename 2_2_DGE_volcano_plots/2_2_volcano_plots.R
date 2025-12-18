#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(grid)   # unit()
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

padj_cut   <- 0.05
lfc_cut    <- 1
max_labels <- 30

# --- VISUAL: NO recortes Y a 80. Solo recortamos X para que se vea el “volcano”,
# --- manteniendo el rango completo de mlog10_padj (~0–310).
xlim_zoom <- c(-7, 7)

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

prep_volcano_df <- function(df) {
  
  needed <- c("log2FoldChange", "padj_safe", "mlog10_padj")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) stop("Faltan columnas requeridas: ", paste(missing, collapse = ", "))
  
  # símbolo de gen: preferir Symbol; si no existe, usar Symbol_autho
  if (!("Symbol" %in% colnames(df))) {
    if ("Symbol_autho" %in% colnames(df)) {
      df <- df %>% mutate(Symbol = as.character(Symbol_autho))
    } else {
      stop("No encuentro columna Symbol ni Symbol_autho.")
    }
  }
  
  df %>%
    mutate(
      Gene = as.character(Symbol),
      log2FoldChange = as.numeric(log2FoldChange),
      padj_used      = as.numeric(padj_safe),
      neglog10_padj  = as.numeric(mlog10_padj)  # ya viene precomputado
    ) %>%
    filter(
      !is.na(Gene), Gene != "",
      is.finite(log2FoldChange),
      !is.na(padj_used), is.finite(padj_used), padj_used > 0,
      !is.na(neglog10_padj), is.finite(neglog10_padj), neglog10_padj >= 0
    ) %>%
    mutate(
      Type = case_when(
        log2FoldChange >=  lfc_cut & padj_used < padj_cut ~ "Upregulated",
        log2FoldChange <= -lfc_cut & padj_used < padj_cut ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      is_adenosine    = Gene %in% adenosine_genes,
      is_adenosine_DE = is_adenosine & padj_used < padj_cut & abs(log2FoldChange) >= lfc_cut
    )
}

make_volcano <- function(res_df, title, subtitle, out_pdf) {
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = Type), size = 0.75, alpha = 0.75) +
    scale_color_manual(values = c(
      "Upregulated"     = "red",
      "Downregulated"   = "blue",
      "Not Significant" = "gray70"
    )) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut),
               color = "black", linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cut),
               color = "black", linetype = "dashed", linewidth = 0.4) +
    labs(
      title    = title,
      subtitle = subtitle,
      x = "log2FoldChange",
      y = expression(-log[10]("padj_safe"))  # (se grafica mlog10_padj)
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank()) +
    # ---- CORRECCIÓN VISUAL ----
  # Solo limitamos X. NO limitamos Y (para que llegue a ~300).
  coord_cartesian(xlim = xlim_zoom)
  
  ad_df <- res_df %>%
    filter(is_adenosine_DE) %>%
    arrange(desc(neglog10_padj)) %>%
    slice_head(n = max_labels)
  
  if (nrow(ad_df) > 0) {
    p <- p +
      geom_point(
        data = ad_df,
        aes(x = log2FoldChange, y = neglog10_padj),
        inherit.aes = FALSE,
        shape  = 21,
        fill   = "purple3",
        color  = "black",
        stroke = 0.35,
        size   = 2.8,
        alpha  = 0.99
      ) +
      ggrepel::geom_label_repel(
        data = ad_df,
        aes(x = log2FoldChange, y = neglog10_padj, label = Gene),
        inherit.aes = FALSE,
        size = 3.4,
        label.size = 0.25,
        label.padding = grid::unit(0.15, "lines"),
        box.padding = 0.40,
        point.padding = 0.22,
        min.segment.length = 0,
        max.overlaps = Inf
      )
  } else {
    message("Nota: No hay genes de adenosina DE (|log2FC|>= ", lfc_cut,
            " y padj_safe<", padj_cut, ") para: ", subtitle)
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

# ===================== SAVE ADENOSINE TABLES ===================== #
# Nota: tu error anterior de select() pasa cuando dplyr no está “ganando” el namespace
# o si 'select' fue enmascarado. Aquí lo dejamos 100% explícito: dplyr::select.

raw_ad <- raw_df %>%
  filter(is_adenosine) %>%
  dplyr::select(Gene, log2FoldChange, padj_used, neglog10_padj, Type, is_adenosine_DE) %>%
  arrange(desc(is_adenosine_DE), Gene)

ae_ad <- ae_df %>%
  filter(is_adenosine) %>%
  dplyr::select(Gene, log2FoldChange, padj_used, neglog10_padj, Type, is_adenosine_DE) %>%
  arrange(desc(is_adenosine_DE), Gene)

readr::write_tsv(raw_ad, file.path(out_dir, "adenosine_genes_status_Ovary_control_GTEx.tsv"))
readr::write_tsv(ae_ad,  file.path(out_dir, "adenosine_genes_status_Reference_control.tsv"))

message("\nDone. Outputs en: ", out_dir)
message("Zoom usado: xlim=", paste(xlim_zoom, collapse = ","),
        "  (Y se deja completo para alcanzar ~300)")
