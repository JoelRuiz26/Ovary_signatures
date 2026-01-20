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
  raw = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  # GEO (estructura que vienes usando: lista con 6 datasets)
  geo = file.path(base_dir, "2_DEG_GEO",  "DEGs_alldsets_GEO.rds")
)

out_dir <- file.path(base_dir, "3_Consensus_DGE_analysis/3_2_DGE_volcano_plots")
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

# ===== NUEVO: GEO (limma) =====
prep_df_geo <- function(df) {
  
  # Gene symbol column (limma outputs típicos + robustez)
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
  
  # logFC column
  lfc_col <- if ("logFC" %in% colnames(df)) "logFC" else if ("log2FoldChange" %in% colnames(df)) "log2FoldChange" else NA
  if (is.na(lfc_col)) stop("GEO: no encuentro columna de logFC (logFC/log2FoldChange).")
  
  # adjusted p-value column
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
      neglog10 = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10)
  
  # En microarreglos puede haber múltiples probes por gen: deja 1 fila por Gene
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
      aes(x = log2FoldChange, y = neglog10, label = Gene),
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
  
  # PDF (vectorial)
  ggsave(out_pdf, p, width = 8, height = 6, device = cairo_pdf)
  
  # PNG 600 dpi
  out_png <- sub("\\.pdf$", ".png", out_pdf)
  ggsave(out_png, p, width = 8, height = 6, dpi = 600, device = "png")
}

# ===================== RUN ===================== #
raw_df <- prep_df(read_rds_safe(paths$raw))
ae_df  <- prep_df(read_rds_safe(paths$ae))

make_volcano(
  raw_df,
  "Ovary tumors vs ovary control GTEx",
  file.path(out_dir, "3_2_1_volcano_Ovary_control_GTEx.pdf")
)

make_volcano(
  ae_df,
  "Ovary tumors vs reference control",
  file.path(out_dir, "3_2_1_volcano_Reference_control_AE.pdf")
)

# ===== NUEVO: GEO volcanos (6 datasets) =====
if (file.exists(paths$geo)) {
  GEO_full <- read_rds_safe(paths$geo)
  
  # si por alguna razón fuera un data.frame, lo convertimos a lista
  if (is.data.frame(GEO_full)) {
    GEO_full <- list(GEO = GEO_full)
  }
  
  if (!is.list(GEO_full) || length(GEO_full) == 0) {
    stop("El objeto GEO no es una lista no-vacía. Revisa: ", paths$geo)
  }
  
  # nombres (idealmente: GSE14407, GSE18520, ...)
  geo_names <- names(GEO_full)
  if (is.null(geo_names) || any(geo_names == "")) {
    geo_names <- paste0("GEO_", seq_along(GEO_full))
  }
  
  for (nm in geo_names) {
    message("GEO volcano: ", nm)
    
    df_geo <- prep_df_geo(GEO_full[[nm]])
    
    # etiqueta limpia para filename
    nm_safe <- gsub("[^A-Za-z0-9_\\-]", "_", nm)
    
    make_volcano(
      df_geo,
      paste0("Ovary tumors vs control (", nm, ")"),
      file.path(out_dir, paste0("3_2_1_volcano_GEO_", nm_safe, ".pdf"))
    )
  }
} else {
  message("Nota: no existe GEO RDS en: ", paths$geo, " (se omiten volcanos GEO).")
}

message("Done. Outputs en: ", out_dir)
