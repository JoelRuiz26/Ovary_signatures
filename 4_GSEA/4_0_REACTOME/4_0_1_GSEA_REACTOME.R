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
  raw = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds")
)

# output folder (Reactome)
out_dir <- file.path(base_dir, "4_GSEA", "4_0_REACTOME")

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

DE_raw <- read_rds_safe(paths$raw)
DE_ae  <- read_rds_safe(paths$ae)

# ===================== REACTOME TERM2GENE ===================== #
# Reactome (MSigDB): C2:CP:REACTOME
react <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  transmute(term = gs_name, entrez = entrez_gene) %>%
  filter(!is.na(entrez)) %>%
  distinct()

TERM2GENE_REACT <- react

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

run_gsea_reactome <- function(geneList, p_cut = 0.01) {
  gsea <- GSEA(
    geneList     = geneList,
    TERM2GENE    = TERM2GENE_REACT,
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

plot_top20_up_down <- function(gsea_df, main_title, subtitle, out_pdf, p_cut = 0.01) {
  if (nrow(gsea_df) == 0) {
    message("No hay resultados para plot: ", subtitle)
    return(invisible(NULL))
  }
  
  dfp <- gsea_df %>%
    filter(!is.na(NES), !is.na(p.adjust), p.adjust > 0) %>%
    mutate(neglog10_fdr = -log10(p.adjust)) %>%
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
    mutate(
      Description = gsub("^REACTOME_", "", Description)
    ) %>%
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
      legend.box = "vertical",
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(face = "bold")
    )
  
  # ====== Guardado PDF + PNG 600 DPI (igual que KEGG) ====== #
  ggsave(out_pdf, p,
         width = 7, height = 8.5, units = "in",
         device = cairo_pdf)
  
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggsave(out_png, p,
         width = 7, height = 8.5, units = "in",
         dpi = 600)
  
  invisible(p)
}

save_outputs <- function(prefix, subtitle, res_list) {
  
  write_tsv(res_list$full,
            file.path(out_dir, paste0(prefix, "_REACTOME_FULL.tsv")))
  write_tsv(res_list$sig,
            file.path(out_dir, paste0(prefix, "_REACTOME_SIG_0.01.tsv")))
  
  plot_top20_up_down(
    gsea_df    = res_list$full,
    main_title = "GSEA Reactome pathways",
    subtitle   = subtitle,
    out_pdf    = file.path(out_dir, paste0(prefix, "_REACTOME_top20.pdf")),
    p_cut      = 0.05
  )
}

# ===================== RUN RAW ===================== #
res_raw <- run_gsea_reactome(make_rank_stat(DE_raw, stat_col = "stat", id_col = "GeneID"))
save_outputs("GTEx", "Ovary tumors vs ovary control GTEx", res_raw)

# ===================== RUN REFERENCE ===================== #
res_ae <- run_gsea_reactome(make_rank_stat(DE_ae, stat_col = "stat", id_col = "GeneID"))
save_outputs("AE", "Ovary tumors vs reference control", res_ae)

message("\nREACTOME GSEA DONE")
