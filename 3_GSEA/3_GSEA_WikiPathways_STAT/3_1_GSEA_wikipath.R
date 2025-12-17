#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw = file.path(base_dir, "0_DGE_raw",         "DE_full_OVARY_DESeq2.rds"),
  ae  = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds")
)

out_dir <- file.path(base_dir, "3_GSEA","3_GSEA_WikiPathways_STAT")

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

DE_raw <- read_rds_safe(paths$raw)
DE_ae  <- read_rds_safe(paths$ae)

# ===================== WIKIPATHWAYS TERM2GENE ===================== #
wp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%
  transmute(term = gs_name, entrez = entrez_gene) %>%
  filter(!is.na(entrez)) %>%
  distinct()

TERM2GENE_WP <- wp

# ===================== HELPERS ===================== #
make_rank_stat <- function(df, stat_col = "stat", id_col = "GeneID") {
  if (!all(c(stat_col, id_col) %in% colnames(df))) {
    stop("Faltan columnas requeridas: ", stat_col, " y/o ", id_col)
  }
  
  r <- df %>%
    transmute(
      entrez = as.character(.data[[id_col]]),
      stat   = as.numeric(.data[[stat_col]])
    ) %>%
    filter(!is.na(entrez), entrez != "", is.finite(stat)) %>%
    group_by(entrez) %>%
    summarise(stat = stat[which.max(abs(stat))], .groups = "drop") %>%
    arrange(desc(stat))
  
  geneList <- r$stat
  names(geneList) <- r$entrez
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
    arrange(NES) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot(top, aes(
    x = NES,
    y = Description,
    color = neglog10_fdr,
    size  = abs(NES)
  )) +
    geom_point(alpha = 0.95, shape = 16) +   # <-- TODO BOLAS (sin triángulos)
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
      plot.subtitle = element_text(size = 10)
    )
  
  ggsave(out_pdf, p, width = 11, height = 8.5, device = cairo_pdf)
  invisible(p)
}

save_outputs <- function(prefix, subtitle, res_list) {
  full_tsv <- file.path(out_dir, paste0(prefix, "_WikiPathways_STAT_FULL.tsv"))
  sig_tsv  <- file.path(out_dir, paste0(prefix, "_WikiPathways_STAT_SIG_FDR0.05.tsv"))
  pdf_top  <- file.path(out_dir, paste0(prefix, "_WikiPathways_STAT_top20_up_down.pdf"))
  
  write_tsv(res_list$full, full_tsv)
  write_tsv(res_list$sig,  sig_tsv)
  
  plot_top20_up_down(
    gsea_df    = res_list$full,
    main_title = "Gene set enrichment analysis using WikiPathways",
    subtitle   = subtitle,
    out_pdf    = pdf_top,
    p_cut      = 0.05
  )
  
  message("OK: ", prefix)
  message("  FULL: ", full_tsv)
  message("  SIG : ", sig_tsv)
  message("  PDF : ", pdf_top)
}

# ===================== RUN (RAW) ===================== #
geneList_raw <- make_rank_stat(DE_raw, stat_col = "stat", id_col = "GeneID")
res_raw <- run_gsea_wp(geneList_raw, p_cut = 0.05)
save_outputs(
  prefix   = "Ovary_control_GTEx",
  subtitle = "Ovary tumors vs ovary control GTEx",
  res_list = res_raw
)

# ===================== RUN (REFERENCE CONTROL / AUTOENCODER) ===================== #
geneList_ae <- make_rank_stat(DE_ae, stat_col = "stat", id_col = "GeneID")
res_ae <- run_gsea_wp(geneList_ae, p_cut = 0.05)
save_outputs(
  prefix   = "Reference_control",
  subtitle = "Ovary tumors vs reference control",
  res_list = res_ae
)

message("\nDone. Outputs en: ", out_dir)
