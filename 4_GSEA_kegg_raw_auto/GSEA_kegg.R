#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
})

# ===================== PATHS ===================== #
base_dir <- "/STORAGE/csbig/jruiz/Ovary_data"

paths <- list(
  raw = file.path(base_dir, "0_DGE_raw", "DE_full_OVARY_DESeq2.rds"),
  ae  = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds")
)

out_dirs <- c(
  "/STORAGE/csbig/jruiz/Ovary_data/3_GSEA_kegg_raw_autoencoder",
  "/home/jruiz/Ovary_signatures/4_GSEA_kegg_raw_auto"
)
invisible(lapply(out_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ===================== SAFE READ ===================== #
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

message("Cargando DESeq2 tables...")
DE_full_raw_DESeq2 <- read_rds_safe(paths$raw)
DE_full_ae_DESeq2  <- read_rds_safe(paths$ae)

# ===================== HELPERS ===================== #
make_rank <- function(df, label = "dataset") {
  
  x <- df %>%
    filter(!is.na(GeneID), !is.na(log2FoldChange)) %>%
    mutate(GeneID = as.character(GeneID)) %>%
    group_by(GeneID) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- x$log2FoldChange
  names(ranks) <- x$GeneID
  ranks <- sort(ranks, decreasing = TRUE)
  
  if (length(ranks) < 1000) {
    warning("Ranking pequeño (", length(ranks), ") en ", label)
  }
  ranks
}

run_gse_kegg <- function(ranks, label,
                         organism = "hsa",
                         pvalueCutoff = 0.05,
                         minGSSize = 10,
                         maxGSSize = 500,
                         eps = 1e-10) {
  
  message("Running gseKEGG for: ", label)
  
  tryCatch({
    gseKEGG(
      geneList     = ranks,
      organism     = organism,
      pvalueCutoff = pvalueCutoff,
      minGSSize    = minGSSize,
      maxGSSize    = maxGSSize,
      eps          = eps,
      verbose      = FALSE
    )
  }, error = function(e) {
    message("ERROR en gseKEGG (", label, "): ", conditionMessage(e))
    NULL
  })
}

save_all <- function(gsea_obj, label) {
  
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("Sin resultados GSEA para: ", label)
    return(invisible(NULL))
  }
  
  df <- as.data.frame(gsea_obj)
  
  # ---------- guardar tablas ----------
  for (od in out_dirs) {
    write_tsv(df, file.path(od, paste0("GSEA_KEGG_", label, ".tsv")))
    saveRDS(gsea_obj, file.path(od, paste0("GSEA_KEGG_", label, ".rds")))
  }
  
  # ---------- plots ----------
  p_dot <- dotplot(gsea_obj, showCategory = 20, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = paste0("GSEA KEGG — ", label))
  
  p_ridge <- ridgeplot(gsea_obj, showCategory = 20) +
    labs(title = paste0("GSEA KEGG (ridge) — ", label))
  
  top_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 20) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  p_nes <- ggplot(top_df, aes(x = Description, y = NES)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top KEGG pathways — ", label),
         x = NULL, y = "NES")
  
  for (od in out_dirs) {
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_dotplot.pdf")),
           p_dot, width = 10, height = 6)
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_ridgeplot.pdf")),
           p_ridge, width = 10, height = 7)
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_NES_top20.pdf")),
           p_nes, width = 10, height = 7)
  }
  
  invisible(df)
}

# ===================== RUN ===================== #
ranks_raw <- make_rank(DE_full_raw_DESeq2, "raw")
ranks_ae  <- make_rank(DE_full_ae_DESeq2,  "autoencoder")

gsea_raw <- run_gse_kegg(ranks_raw, "raw")
gsea_ae  <- run_gse_kegg(ranks_ae,  "autoencoder")

df_raw <- save_all(gsea_raw, "raw")
df_ae  <- save_all(gsea_ae,  "autoencoder")

message("DONE.")
