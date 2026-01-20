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

# Meta-volcano settings
min_consistent_n <- 7  # genes presentes en >=5/8 datasets para evaluar consistencia
p_floor <- 1e-300      # evitar log(0) en Fisher y -log10

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
  
  ggsave(out_pdf, p, width = 8, height = 6, device = cairo_pdf)
  out_png <- sub("\\.pdf$", ".png", out_pdf)
  ggsave(out_png, p, width = 8, height = 6, dpi = 600, device = "png")
}

# ===================== META HELPERS ===================== #
combine_p_fisher <- function(p_vec) {
  p_vec <- p_vec[is.finite(p_vec) & !is.na(p_vec)]
  if (length(p_vec) == 0) return(NA_real_)
  p_vec <- pmax(p_vec, p_floor)
  stat <- -2 * sum(log(p_vec))
  df <- 2 * length(p_vec)
  pchisq(stat, df = df, lower.tail = FALSE)
}

# Construye tabla long con Dataset/Gene/log2FC/padj/neglog10
build_all_long <- function(raw_df, ae_df, GEO_full) {
  out <- list()
  
  out[["GTEx"]] <- raw_df %>%
    transmute(Dataset = "GTEx", Gene, log2FoldChange, padj_used, neglog10)
  
  out[["AE"]] <- ae_df %>%
    transmute(Dataset = "AE", Gene, log2FoldChange, padj_used, neglog10)
  
  if (!is.null(GEO_full) && is.list(GEO_full) && length(GEO_full) > 0) {
    geo_names <- names(GEO_full)
    if (is.null(geo_names) || any(geo_names == "")) {
      geo_names <- paste0("GEO_", seq_along(GEO_full))
    }
    for (nm in geo_names) {
      df_geo <- prep_df_geo(GEO_full[[nm]])
      out[[paste0("GEO_", nm)]] <- df_geo %>%
        transmute(Dataset = paste0("GEO_", nm), Gene, log2FoldChange, padj_used, neglog10)
    }
  }
  
  bind_rows(out)
}

make_meta_volcano_consistent <- function(all_long_df, out_pdf) {
  
  # 1) RANGOS por gen en todos los datasets + consistencia de dirección
  meta_ranges <- all_long_df %>%
    group_by(Gene) %>%
    summarise(
      n_datasets = n_distinct(Dataset),
      fc_min = min(log2FoldChange, na.rm = TRUE),
      fc_max = max(log2FoldChange, na.rm = TRUE),
      nl_min = min(neglog10, na.rm = TRUE),
      nl_max = max(neglog10, na.rm = TRUE),
      all_pos = all(log2FoldChange > 0, na.rm = TRUE),
      all_neg = all(log2FoldChange < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dir = case_when(
        all_pos ~ "Up",
        all_neg ~ "Down",
        TRUE ~ "Mixed"
      ),
      common_consistent = (n_datasets >= min_consistent_n) & (dir %in% c("Up", "Down")),
      color_class = case_when(
        common_consistent & dir == "Up"   ~ "Common Up",
        common_consistent & dir == "Down" ~ "Common Down",
        TRUE ~ "Other"
      ),
      is_adenosine = Gene %in% adenosine_genes
    )
  
  # 2) PUNTO ÚNICO POR GEN (evitar duplicados):
  #    - x: logFC con mayor |logFC|
  #    - y: neglog10 máximo (mayor significancia)
  #    Nota: se calcula desde all_long_df (punto "más extremo" observado)
  meta_point <- all_long_df %>%
    group_by(Gene) %>%
    summarise(
      x_fc = log2FoldChange[which.max(abs(log2FoldChange))],
      y_nl = max(neglog10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_nl = pmin(y_nl, max_neglog10)
    )
  
  meta <- meta_ranges %>%
    inner_join(meta_point, by = "Gene")
  
  meta_hi <- meta %>% filter(common_consistent)
  
  # --- recorte visual para que no haya "p-values por las nubes" ---
  y_cap <- 60  # ajusta si quieres (40-80 suele verse bien)
  
  meta_plot <- meta %>%
    mutate(y_plot = pmin(y_nl, y_cap))
  
  meta_hi <- meta_plot %>%
    filter(common_consistent, color_class != "Other")
  
  p <- ggplot(meta_plot, aes(x = x_fc, y = y_plot)) +
    # Fondo (todos) como volcano tradicional
    geom_point(color = "gray80", size = 0.35, alpha = 0.25) +
    
    # Resaltados (comunes-consistentes) arriba
    geom_point(
      data = meta_hi,
      aes(color = color_class),
      size = 1.15,
      alpha = 0.95
    ) +
    
    # Umbrales (referencia)
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.35) +
    geom_hline(yintercept = -log10(padj_cut),     linetype = "dashed", linewidth = 0.35) +
    
    # Etiquetas: solo adenosina dentro de los resaltados
    ggrepel::geom_label_repel(
      data = meta_hi %>% filter(is_adenosine),
      aes(label = Gene),
      size = 3.1,
      max.overlaps = Inf,
      fill = "white",
      color = "black",
      label.size = 0.25,
      min.segment.length = 0
    ) +
    
    scale_color_manual(values = c(
      "Common Up"   = "red",
      "Common Down" = "blue"
    )) +
    
    scale_y_continuous(
      limits = c(0, y_cap),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    
    labs(
      title = "Meta-volcano (GTEx + AE + GEO)",
      subtitle = paste0(
        "Gray: all genes (1 point/gene; y capped at ", y_cap, "). ",
        "Colored: genes present in ≥", min_consistent_n, " datasets with consistent direction."
      ),
      x = "log2FC (picked as max |log2FC| across datasets)",
      y = "-log10(padj) (picked as max across datasets; capped for display)",
      color = "Common genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  
  ggsave(out_pdf, p, width = 9.2, height = 6.6, device = cairo_pdf)
  out_png <- sub("\\.pdf$", ".png", out_pdf)
  ggsave(out_png, p, width = 9.2, height = 6.6, dpi = 600, device = "png")
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

# ===== GEO volcanos (6 datasets) =====
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
    
    make_volcano(
      df_geo,
      paste0("Ovary tumors vs control (", nm, ")"),
      file.path(out_dir, paste0("3_2_1_volcano_GEO_", nm_safe, ".pdf"))
    )
  }
} else {
  message("Nota: no existe GEO RDS en: ", paths$geo, " (se omiten volcanos GEO).")
}

# ===== NUEVO: META-VOLCANO CONSISTENTE =====
all_long_df <- build_all_long(raw_df, ae_df, GEO_full)

make_meta_volcano_consistent(
  all_long_df,
  file.path(out_dir, "3_2_2_META_volcano_consistent_GTEx_AE_GEO.pdf")
)

message("Done. Outputs en: ", out_dir)
