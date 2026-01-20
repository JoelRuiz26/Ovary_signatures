#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

# Input: lista con 6 datasets (GSE...)
geo_path <- file.path(base_dir, "2_DEG_GEO", "DEGs_alldsets_GEO.rds")

# Output dir (MISMO que tu script)
out_dir <- file.path(base_dir, "4_GSEA", "4_0_KEGG")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Threshold (cambio pedido)
padj_cut <- 0.05

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

GEO_full <- read_rds_safe(geo_path)
if (!is.list(GEO_full) || length(GEO_full) == 0) {
  stop("El objeto GEO_full no es una lista válida o está vacía.")
}

# ===================== HELPERS ===================== #

# Map SYMBOL -> ENTREZID (necesario para gseKEGG)
symbols_to_entrez <- function(symbols_vec) {
  symbols_vec <- unique(na.omit(as.character(symbols_vec)))
  symbols_vec <- symbols_vec[symbols_vec != ""]
  if (length(symbols_vec) == 0) return(tibble(SYMBOL = character(), ENTREZID = character()))
  
  suppressWarnings({
    mm <- bitr(
      symbols_vec,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  })
  
  if (is.null(mm) || nrow(mm) == 0) {
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  }
  
  mm %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    transmute(
      SYMBOL   = as.character(SYMBOL),
      ENTREZID = as.character(ENTREZID)
    )
}

# make_rank_stat para GEO (ranking por t) + impresión antes/después de anotar
make_rank_stat_geo <- function(df,
                               symbol_col = "Gene.symbol",
                               stat_col   = "t",
                               verbose    = FALSE) {
  
  if (!all(c(symbol_col, stat_col) %in% colnames(df))) {
    stop("Faltan columnas requeridas: ", symbol_col, " y/o ", stat_col)
  }
  
  tmp <- df %>%
    transmute(
      SYMBOL = as.character(.data[[symbol_col]]),
      stat   = as.numeric(.data[[stat_col]])
    ) %>%
    filter(!is.na(SYMBOL), SYMBOL != "", is.finite(stat))
  
  n_sym_unique <- dplyr::n_distinct(tmp$SYMBOL)
  
  mm <- symbols_to_entrez(tmp$SYMBOL)
  if (nrow(mm) == 0) stop("No se pudo mapear ningún SYMBOL a ENTREZID.")
  
  n_sym_mapped <- dplyr::n_distinct(mm$SYMBOL)
  map_pct <- if (n_sym_unique == 0) 0 else round(100 * n_sym_mapped / n_sym_unique, 2)
  
  r <- tmp %>%
    inner_join(mm, by = "SYMBOL") %>%
    transmute(
      entrez = as.character(ENTREZID),
      stat   = as.numeric(stat)
    ) %>%
    filter(!is.na(entrez), entrez != "", is.finite(stat)) %>%
    group_by(entrez) %>%
    summarise(stat = stat[which.max(abs(stat))], .groups = "drop") %>%
    arrange(desc(stat))
  
  geneList <- r$stat
  names(geneList) <- r$entrez
  
  if (verbose) {
    message("  Genes (rows) en tabla: ", nrow(df))
    message("  SYMBOL únicos (pre-mapeo): ", n_sym_unique)
    message("  SYMBOL que mapearon a ENTREZ: ", n_sym_mapped, " (", map_pct, "%)")
    message("  ENTREZ únicos que entran a GSEA: ", length(geneList))
    message("  Ejemplo ENTREZ (primeros 10): ", paste(head(names(geneList), 10), collapse = ", "))
  }
  
  geneList
}

run_gsea_kegg <- function(geneList, p_cut = 0.05) {
  gsea <- gseKEGG(
    geneList     = geneList,
    organism     = "hsa",
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  
  gsea_df  <- as.data.frame(gsea)
  gsea_sig <- gsea_df %>% filter(!is.na(p.adjust), p.adjust < p_cut)
  
  list(gsea = gsea, full = gsea_df, sig = gsea_sig)
}

plot_top20_up_down <- function(gsea_df, main_title, subtitle, out_pdf, p_cut = 0.05) {
  
  dfp <- gsea_df %>%
    filter(!is.na(NES), !is.na(p.adjust), p.adjust < p_cut) %>%
    mutate(neglog10_fdr = -log10(p.adjust))
  
  top_up <- dfp %>%
    filter(NES > 0) %>%
    arrange(p.adjust, desc(NES)) %>%
    slice_head(n = 20)
  
  top_down <- dfp %>%
    filter(NES < 0) %>%
    arrange(p.adjust, NES) %>%
    slice_head(n = 20)
  
  top <- bind_rows(top_up, top_down) %>%
    arrange(NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot(top, aes(
    x = NES,
    y = Description,
    color = neglog10_fdr,
    size  = abs(NES)
  )) +
    geom_point(shape = 16, alpha = 0.95) +
    scale_color_viridis_c(name = expression(-log[10]("FDR"))) +
    scale_size_continuous(name = "|NES|") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = main_title,
      subtitle = subtitle,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "bold")
    )
  
  ggsave(out_pdf, p,
         width = 7, height = 8.5, units = "in",
         device = cairo_pdf)
  
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggsave(out_png, p,
         width = 7, height = 8.5, units = "in",
         dpi = 600)
}

save_outputs <- function(prefix, subtitle, res_list, p_cut = 0.05) {
  
  write_tsv(res_list$full,
            file.path(out_dir, paste0(prefix, "_KEGG_FULL.tsv")))
  write_tsv(res_list$sig,
            file.path(out_dir, paste0(prefix, "_KEGG_SIG_0.05.tsv")))
  
  plot_top20_up_down(
    gsea_df    = res_list$full,
    main_title = "GSEA KEGG pathways",
    subtitle   = subtitle,
    out_pdf    = file.path(out_dir, paste0(prefix, "_KEGG_top20.pdf")),
    p_cut      = p_cut
  )
}

# ===================== RUN GEO (6 datasets) ===================== #
for (dset_name in names(GEO_full)) {
  
  message("\n==============================")
  message("Running KEGG GSEA for: ", dset_name)
  message("==============================")
  
  df <- GEO_full[[dset_name]]
  if (!is.data.frame(df)) {
    warning("Saltando ", dset_name, ": no es data.frame")
    next
  }
  
  # Ranking por t; imprime antes/después de anotar
  geneList <- make_rank_stat_geo(df, symbol_col = "Gene.symbol", stat_col = "t", verbose = TRUE)
  
  res <- run_gsea_kegg(geneList, p_cut = padj_cut)
  
  message("  KEGG total rows: ", nrow(res$full))
  message("  KEGG sig (FDR<", padj_cut, "): ", nrow(res$sig))
  
  prefix   <- paste0("GEO_", dset_name)
  subtitle <- paste0("Ovary tumors vs control (", dset_name, ")")
  
  save_outputs(prefix, subtitle, res, p_cut = padj_cut)
}

message("\nKEGG GSEA DONE (GEO 6 datasets)")
