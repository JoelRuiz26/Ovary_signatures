#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures/4_GSEA"

dbs <- list(
  KEGG         = "4_0_KEGG",
  Reactome     = "4_0_REACTOME",
  WikiPathways = "4_0_WIKIPATHWAYS"
)

fdr_tag <- "0.05"      # tus SIG_0.05.tsv
min_k   <- 5           # al menos 5 datasets (de 8)
top_n   <- 20          # top 20 términos

# ===================== HELPERS ===================== #
load_sig <- function(dir, prefix, fdr_tag = "0.05") {
  
  files <- list.files(
    path = dir,
    pattern = paste0("^", prefix, ".*_SIG_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    files <- list.files(
      path = dir,
      pattern = paste0("^", prefix, ".*_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
      full.names = TRUE
    )
  }
  
  if (length(files) == 0) {
    stop("No se encontró archivo SIG para: ", prefix, " (FDR ", fdr_tag, ") en ", dir)
  }
  if (length(files) > 1) {
    message("Más de un archivo encontrado para ", prefix, " en ", dir, "\nUsando: ", files[1])
  }
  
  readr::read_tsv(files[1], show_col_types = FALSE)
}

clean_description <- function(df, db_name) {
  if (db_name == "Reactome") {
    df <- df %>% mutate(Description = gsub("^REACTOME_", "", Description))
  }
  if (db_name == "WikiPathways") {
    df <- df %>% mutate(Description = gsub("^WP_", "", Description))
  }
  df
}

# ---- UpSet plot (intenta ComplexUpset; si no, usa UpSetR) ----
plot_upset_safe <- function(inc_df, set_cols, title_txt, out_pdf, out_png) {
  
  # inc_df: data.frame con columnas lógicas/0-1 por set, y (opcional) otras cols
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    # ComplexUpset usa ggplot2
    p <- ComplexUpset::upset(
      inc_df,
      intersect = set_cols,
      name = "Intersection",
      width_ratio = 0.15
    ) +
      ggplot2::labs(title = title_txt)
    
    ggplot2::ggsave(out_pdf, p, width = 12, height = 6, units = "in", device = cairo_pdf)
    ggplot2::ggsave(out_png, p, width = 12, height = 6, units = "in", dpi = 600)
    return(invisible(TRUE))
  }
  
  if (requireNamespace("UpSetR", quietly = TRUE)) {
    grDevices::pdf(out_pdf, width = 12, height = 6)
    UpSetR::upset(
      inc_df,
      sets = set_cols,
      keep.order = TRUE,
      order.by = "freq",
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size"
    )
    grDevices::dev.off()
    
    # PNG (base)
    grDevices::png(out_png, width = 12, height = 6, units = "in", res = 600)
    UpSetR::upset(
      inc_df,
      sets = set_cols,
      keep.order = TRUE,
      order.by = "freq",
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size"
    )
    grDevices::dev.off()
    return(invisible(TRUE))
  }
  
  warning("No está instalado ComplexUpset ni UpSetR. Instala uno:\n",
          "  install.packages('ComplexUpset')  # recomendado\n",
          "  # o\n  install.packages('UpSetR')\n",
          "No se generó el UpSet plot.")
  invisible(FALSE)
}

# Bubble plot para top shared terms
plot_bubble_top <- function(top_tbl, title_txt, subtitle_txt, out_pdf, out_png) {
  
  if (nrow(top_tbl) == 0) {
    message("No hay términos para bubble plot: ", title_txt)
    return(invisible(NULL))
  }
  
  dfp <- top_tbl %>%
    mutate(
      neglog10_fdr = -log10(padj_min),
      Description = factor(Description, levels = rev(Description))
    )
  
  p <- ggplot(dfp, aes(
    x = NES_mean,
    y = Description,
    color = neglog10_fdr,
    size = n_datasets
  )) +
    geom_point(alpha = 0.95, shape = 16) +
    scale_color_viridis_c(name = expression(-log[10]("min FDR"))) +
    scale_size_continuous(name = "N datasets") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Mean NES across datasets (present only)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold")
    )
  
  ggsave(out_pdf, p, width = 11, height = 8.5, units = "in", device = cairo_pdf)
  ggsave(out_png, p, width = 11, height = 8.5, units = "in", dpi = 600)
  
  invisible(p)
}

# ===================== MAIN LOOP ===================== #
for (db in names(dbs)) {
  
  message("\n==============================")
  message("Processing ", db, " (UpSet + top shared >= ", min_k, ") ...")
  message("==============================")
  
  db_dir <- file.path(base_dir, dbs[[db]])
  shared_dir <- file.path(db_dir, "shared")
  dir.create(shared_dir, showWarnings = FALSE)
  
  # --------- Detectar GEO prefijos (robusto) ---------
  geo_sig_files <- list.files(
    path = db_dir,
    pattern = paste0("^GEO_.*_SIG_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
    full.names = FALSE
  )
  
  # Extrae "GEO_GSE####" de: GEO_GSE14407_REACTOME_SIG_0.05.tsv (etc.)
  geo_prefixes <- character()
  if (length(geo_sig_files) > 0) {
    geo_prefixes <- unique(sub(paste0("_[A-Z]+_SIG_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"), "", geo_sig_files))
    # por seguridad, filtra solo los que realmente empiezan con GEO_
    geo_prefixes <- geo_prefixes[grepl("^GEO_", geo_prefixes)]
  }
  
  dataset_prefixes <- c("GTEx", "AE", sort(geo_prefixes))
  message("Datasets detectados: ", paste(dataset_prefixes, collapse = ", "))
  
  # --------- Cargar todos los SIG de los 8 datasets ---------
  all_long <- list()
  
  for (pref in dataset_prefixes) {
    df <- load_sig(db_dir, pref, fdr_tag = fdr_tag) %>%
      clean_description(., db) %>%
      transmute(
        Dataset = pref,
        Description = as.character(Description),
        NES = as.numeric(NES),
        p.adjust = as.numeric(p.adjust)
      ) %>%
      filter(!is.na(Description), Description != "", !is.na(NES), is.finite(NES),
             !is.na(p.adjust), is.finite(p.adjust))
    
    all_long[[pref]] <- df
  }
  
  all_long_df <- bind_rows(all_long)
  
  # --------- Incidencia (término presente / ausente por dataset) ---------
  inc <- all_long_df %>%
    distinct(Dataset, Description) %>%
    mutate(present = 1L) %>%
    tidyr::pivot_wider(
      names_from = Dataset,
      values_from = present,
      values_fill = 0L
    )
  
  # Si tidyr no está, instrucción clara (pero no rompas silenciosamente)
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Falta el paquete 'tidyr' para construir la matriz de incidencia. Instala:\n",
         "  install.packages('tidyr')\n")
  }
  
  # tidyr pivot_wider ya se usó arriba; si no existía, ya habría fallado.
  # (lo dejamos así para que el error sea explícito)
  
  # Calcula cuántos datasets contienen cada término
  set_cols <- dataset_prefixes
  inc <- inc %>%
    mutate(n_datasets = rowSums(across(all_of(set_cols))))
  
  # Filtra términos compartidos por >= min_k datasets
  inc_k <- inc %>% filter(n_datasets >= min_k)
  
  message("Términos con presencia en >= ", min_k, " datasets: ", nrow(inc_k))
  
  # --------- UpSet plot (solo términos >= min_k) ---------
  upset_pdf <- file.path(shared_dir, paste0(db, "_UpSet_shared_ge_", min_k, "_of_", length(set_cols), ".pdf"))
  upset_png <- sub("\\.pdf$", ".png", upset_pdf, ignore.case = TRUE)
  
  # Asegura columnas binarias (0/1) para UpSetR/ComplexUpset
  inc_plot <- inc_k %>%
    select(all_of(set_cols)) %>%
    mutate(across(everything(), ~ as.integer(.x > 0)))
  
  plot_upset_safe(
    inc_df   = inc_plot,
    set_cols = set_cols,
    title_txt = paste0(db, " | terms shared in \u2265", min_k, " / ", length(set_cols), " datasets (FDR<", fdr_tag, ")"),
    out_pdf  = upset_pdf,
    out_png  = upset_png
  )
  
  # --------- Tabla resumen para términos >= min_k y seleccionar Top 20 ---------
  # Resumen por término: frecuencia, min FDR, mean NES (solo donde aparece)
  term_summary <- all_long_df %>%
    semi_join(inc_k %>% select(Description), by = "Description") %>%
    group_by(Description) %>%
    summarise(
      n_datasets = n_distinct(Dataset),
      padj_min   = min(p.adjust, na.rm = TRUE),
      NES_mean   = mean(NES, na.rm = TRUE),
      NES_abs_mean = mean(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_datasets), padj_min, desc(NES_abs_mean))
  
  top_terms <- term_summary %>% slice_head(n = top_n)
  
  out_tsv <- file.path(shared_dir, paste0(db, "_Top", top_n, "_terms_shared_ge_", min_k, "_datasets.tsv"))
  write_tsv(top_terms, out_tsv)
  
  # --------- Bubble plot (solo top20 de >= min_k datasets) ---------
  bub_pdf <- file.path(shared_dir, paste0(db, "_Top", top_n, "_bubble_shared_ge_", min_k, ".pdf"))
  bub_png <- sub("\\.pdf$", ".png", bub_pdf, ignore.case = TRUE)
  
  plot_bubble_top(
    top_tbl = top_terms,
    title_txt = paste0("Top ", top_n, " ", db, " terms shared in \u2265", min_k, " datasets"),
    subtitle_txt = paste0("Ranking: n_datasets desc, min FDR asc, |NES| mean desc (FDR<", fdr_tag, ")"),
    out_pdf = bub_pdf,
    out_png = bub_png
  )
  
  message("Saved:")
  message("  UpSet: ", upset_pdf)
  message("  Top terms TSV: ", out_tsv)
  message("  Bubble: ", bub_pdf)
}

message("\nDONE: UpSet + Top terms + Bubble plots for KEGG/Reactome/WikiPathways.")
