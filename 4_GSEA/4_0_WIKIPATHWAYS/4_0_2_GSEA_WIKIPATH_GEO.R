#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(org.Hs.eg.db)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

# Input: lista con 6 datasets (GSE...)
geo_path <- file.path(base_dir, "2_DEG_GEO", "DEGs_alldsets_GEO.rds")

# Output folder (consistent with 4_GSEA convention)
out_dir <- file.path(base_dir, "4_GSEA", "4_0_WIKIPATHWAYS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Threshold (FDR)
p_cut <- 0.05

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

GEO_full <- read_rds_safe(geo_path)
if (!is.list(GEO_full) || length(GEO_full) == 0) {
  stop("El objeto GEO_full no es una lista válida o está vacía.")
}

# ===================== WIKIPATHWAYS TERM2GENE ===================== #
wp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
  transmute(term = gs_name, entrez = entrez_gene) %>%
  filter(!is.na(entrez)) %>%
  distinct()

TERM2GENE_WP <- wp

# ===================== HELPERS ===================== #

# Map SYMBOL -> ENTREZID
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

run_gsea_wp <- function(geneList, p_cut = 0.05) {
  gsea <- GSEA(
    geneList     = geneList,
    TERM2GENE    = TERM2GENE_WP,
    pvalueCutoff = 1,       # guardar FULL; filtrar después
    minGSSize    = 10,
    maxGSSize    = 500,
    eps          = 0,
    verbose      = FALSE
  )
  
  gsea_df  <- as.data.frame(gsea)
  gsea_sig <- gsea_df %>% filter(!is.na(p.adjust), p.adjust < p_cut)
  
  list(gsea = gsea, full = gsea_df, sig = gsea_sig)
}

plot_top20_up_down <- function(gsea_df, main_title, subtitle, out_pdf, p_cut = 0.05) {
  if (nrow(gsea_df) == 0) {
    message("No hay resultados para plot: ", subtitle)
    return(invisible(NULL))
  }
  
  dfp <- gsea_df %>%
    filter(!is.na(NES), !is.na(p.adjust), p.adjust > 0) %>%
    mutate(
      neglog10_fdr = -log10(p.adjust),
      Direction = ifelse(NES >= 0, "Upregulated pathways", "Downregulated pathways")
    ) %>%
    filter(p.adjust < p_cut)
  
  if (nrow(dfp) == 0) {
    message("No hay vías significativas (FDR < ", p_cut, ") para plot: ", subtitle)
    return(invisible(NULL))
  }
  
  top_up <- dfp %>%
    filter(NES > 0) %>%
    arrange(p.adjust, desc(NES)) %>%
    slice_head(n = 20)
  
  top_down <- dfp %>%
    filter(NES < 0) %>%
    arrange(p.adjust, NES) %>%
    slice_head(n = 20)
  
  top <- bind_rows(top_up, top_down) %>%
    mutate(Description = gsub("^WP_", "", Description)) %>%
    arrange(NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot(top, aes(
    x = NES,
    y = Description,
    color = neglog10_fdr,
    size  = abs(NES)
  )) +
    geom_point(alpha = 0.95, shape = 16) +
    scale_color_viridis_c(name = expression(-log[10]("FDR"))) +
    scale_size_continuous(name = "|NES|") +
    labs(
      title = main_title,
      subtitle = subtitle,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.box = "vertical",
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(face = "bold")
    )
  
  ggsave(out_pdf, p, width = 11, height = 8.5, units = "in", device = cairo_pdf)
  
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggsave(out_png, p, width = 11, height = 8.5, units = "in", dpi = 600)
  
  invisible(p)
}

save_outputs <- function(prefix, subtitle, res_list, p_cut = 0.05) {
  
  full_tsv <- file.path(out_dir, paste0(prefix, "_WIKIPATHWAYS_FULL.tsv"))
  sig_tsv  <- file.path(out_dir, paste0(prefix, "_WIKIPATHWAYS_SIG_0.05.tsv"))
  pdf_top  <- file.path(out_dir, paste0(prefix, "_WIKIPATHWAYS_top20.pdf"))
  
  write_tsv(res_list$full, full_tsv)
  write_tsv(res_list$sig,  sig_tsv)
  
  plot_top20_up_down(
    gsea_df    = res_list$full,
    main_title = "GSEA WikiPathways",
    subtitle   = subtitle,
    out_pdf    = pdf_top,
    p_cut      = p_cut
  )
  
  message("OK: ", prefix)
  message("  FULL: ", full_tsv)
  message("  SIG : ", sig_tsv)
  message("  PDF : ", pdf_top)
  message("  PNG : ", sub("\\.pdf$", ".png", pdf_top, ignore.case = TRUE))
}

# ===================== RUN GEO (6 datasets) ===================== #
for (dset_name in names(GEO_full)) {
  
  message("\n==============================")
  message("Running WIKIPATHWAYS GSEA for: ", dset_name)
  message("==============================")
  
  df <- GEO_full[[dset_name]]
  if (!is.data.frame(df)) {
    warning("Saltando ", dset_name, ": no es data.frame")
    next
  }
  
  geneList <- make_rank_stat_geo(df, symbol_col = "Gene.symbol", stat_col = "t", verbose = TRUE)
  
  res <- run_gsea_wp(geneList, p_cut = p_cut)
  
  message("  WIKIPATHWAYS total rows: ", nrow(res$full))
  message("  WIKIPATHWAYS sig (FDR<", p_cut, "): ", nrow(res$sig))
  
  prefix   <- paste0("GEO_", dset_name)
  subtitle <- paste0("Ovary tumors vs control (", dset_name, ")")
  
  save_outputs(prefix, subtitle, res, p_cut = p_cut)
}

message("\nWIKIPATHWAYS GSEA DONE (GEO 6 datasets). Outputs en: ", out_dir)
