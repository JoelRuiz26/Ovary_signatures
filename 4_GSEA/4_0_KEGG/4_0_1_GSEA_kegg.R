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

paths <- list(
  raw = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds")
)

out_dir <- file.path(base_dir, "4_GSEA", "4_0_KEGG")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

DE_raw <- read_rds_safe(paths$raw)
DE_ae  <- read_rds_safe(paths$ae)

# ===================== HELPERS ===================== #
make_rank_stat <- function(df, stat_col = "stat", id_col = "GeneID") {
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
    mutate(
      neglog10_fdr = -log10(p.adjust)
    )
  
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
  
  # ====== CAMBIOS MÍNIMOS AQUÍ ====== #
  ggsave(out_pdf, p,
         width = 7, height = 8.5, units = "in",
         device = cairo_pdf)
  
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggsave(out_png, p,
         width = 7, height = 8.5, units = "in",
         dpi = 600)
}

save_outputs <- function(prefix, subtitle, res_list) {
  
  write_tsv(res_list$full,
            file.path(out_dir, paste0(prefix, "_KEGG_FULL.tsv")))
  write_tsv(res_list$sig,
            file.path(out_dir, paste0(prefix, "_KEGG_SIG_0.05.tsv")))
  
  plot_top20_up_down(
    gsea_df    = res_list$full,
    main_title = "GSEA KEGG pathways",
    subtitle   = subtitle,
    out_pdf    = file.path(out_dir, paste0(prefix, "_KEGG_top20.pdf"))
  )
}

# ===================== RUN RAW ===================== #
res_raw <- run_gsea_kegg(make_rank_stat(DE_raw))
save_outputs("GTEx", "Ovary tumors vs ovary control GTEx", res_raw)
ls
# ===================== RUN REFERENCE ===================== #
res_ae <- run_gsea_kegg(make_rank_stat(DE_ae))
save_outputs("AE", "Ovary tumors vs reference control", res_ae)

message("\nKEGG GSEA DONE")
