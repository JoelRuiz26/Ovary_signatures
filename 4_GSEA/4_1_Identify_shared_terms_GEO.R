#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
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

# Barplot (ordenado por NES_mean: mayor -> menor)
plot_bar_top <- function(top_tbl, title_txt, subtitle_txt, out_pdf, out_png) {
  
  if (nrow(top_tbl) == 0) {
    message("No hay términos para barplot: ", title_txt)
    return(invisible(NULL))
  }
  
  # ORDEN CORRECTO: NES_mean descendente (visual y estadístico)
  dfp <- top_tbl %>%
    mutate(
      neglog10_fdr = -log10(padj_min)
    ) %>%
    arrange(desc(NES_mean)) %>%
    mutate(
      Description = factor(Description, levels = Description)
    )
  
  p <- ggplot(dfp, aes(
    x = NES_mean,
    y = Description,
    fill = neglog10_fdr
  )) +
    geom_col(width = 0.75, alpha = 0.95) +
    scale_fill_viridis_c(name = expression(-log[10]("min FDR"))) +
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
  
  ggsave(out_pdf, p,
         width = 11, height = 8.5, units = "in",
         device = cairo_pdf)
  
  ggsave(out_png, p,
         width = 11, height = 8.5, units = "in",
         dpi = 600)
  
  invisible(p)
}


# ===================== MAIN LOOP ===================== #
for (db in names(dbs)) {
  
  message("\n==============================")
  message("Processing ", db, " (top shared >= ", min_k, ") ...")
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
  
  geo_prefixes <- character()
  if (length(geo_sig_files) > 0) {
    geo_prefixes <- unique(sub(
      paste0("_[A-Z]+_SIG_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
      "",
      geo_sig_files
    ))
    geo_prefixes <- geo_prefixes[grepl("^GEO_", geo_prefixes)]
  }
  
  # Orden deseado: GEO (ordenado) -> GTEx -> AE
  dataset_prefixes <- c(sort(geo_prefixes), "GTEx", "AE")
  message("Datasets detectados (orden): ", paste(dataset_prefixes, collapse = ", "))
  
  # --------- Cargar todos los SIG ---------
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
      filter(
        !is.na(Description), Description != "",
        !is.na(NES), is.finite(NES),
        !is.na(p.adjust), is.finite(p.adjust)
      )
    
    all_long[[pref]] <- df
  }
  
  all_long_df <- bind_rows(all_long)
  
  # --------- Incidencia para filtrar términos en >= min_k datasets ---------
  inc <- all_long_df %>%
    distinct(Dataset, Description) %>%
    mutate(present = 1L) %>%
    tidyr::pivot_wider(
      names_from = Dataset,
      values_from = present,
      values_fill = 0L
    )
  
  set_cols <- dataset_prefixes
  
  inc <- inc %>%
    mutate(n_datasets = rowSums(across(all_of(set_cols))))
  
  inc_k <- inc %>% filter(n_datasets >= min_k)
  
  message("Términos con presencia en >= ", min_k, " datasets: ", nrow(inc_k))
  
  # --------- Tabla resumen + Top 20 ---------
  term_summary <- all_long_df %>%
    semi_join(inc_k %>% select(Description), by = "Description") %>%
    group_by(Description) %>%
    summarise(
      n_datasets   = n_distinct(Dataset),
      padj_min     = min(p.adjust, na.rm = TRUE),
      NES_mean     = mean(NES, na.rm = TRUE),
      NES_abs_mean = mean(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_datasets), padj_min, desc(NES_abs_mean))
  
  top_terms <- term_summary %>% slice_head(n = top_n)
  
  out_tsv <- file.path(shared_dir, paste0(db, "_Top", top_n, "_terms_shared_ge_", min_k, "_datasets.tsv"))
  write_tsv(top_terms, out_tsv)
  
  # --------- Barplot (ordenado) ---------
  bar_pdf <- file.path(shared_dir, paste0(db, "_Top", top_n, "_bar_shared_ge_", min_k, ".pdf"))
  bar_png <- sub("\\.pdf$", ".png", bar_pdf, ignore.case = TRUE)
  
  plot_bar_top(
    top_tbl = top_terms,
    title_txt = paste0("Top ", top_n, " ", db, " terms shared in ≥", min_k, " datasets"),
    subtitle_txt = paste0("Ranking: n_datasets desc, min FDR asc, |NES| mean desc (FDR<", fdr_tag,
                          "). Bars ordered by NES_mean."),
    out_pdf = bar_pdf,
    out_png = bar_png
  )
  
  message("Saved:")
  message("  Top terms TSV: ", out_tsv)
  message("  Barplot: ", bar_pdf)
}

message("\nDONE: Top terms + Barplots for KEGG/Reactome/WikiPathways (GEO+GTEx+AE).")
