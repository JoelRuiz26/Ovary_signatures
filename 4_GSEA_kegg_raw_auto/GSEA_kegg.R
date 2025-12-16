#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
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

# ===================== PARAMETERS ===================== #
# C6 = Oncogenic signatures (MSigDB)
MSIG_CATEGORY    <- "C6"
MSIG_SUBCATEGORY <- NULL  # no aplica para C6

# GSEA settings
MIN_SIZE <- 15
MAX_SIZE <- 500
NPERM    <- 20000   # más estable, menos warnings tipo rmSimple

# Output settings
SHOW_TOP <- 30      # más categorías en plots (y PDFs más altos)

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

# C6 oncogenic sets from MSigDB (Entrez IDs)
get_msig_oncogenic <- function() {
  message("Cargando MSigDB ", MSIG_CATEGORY, " (oncogenic signatures)...")
  msig <- msigdbr(species = "Homo sapiens", category = MSIG_CATEGORY)
  
  # pathways: named list of Entrez IDs (character)
  pathways <- msig %>%
    transmute(gs_name, entrez_gene) %>%
    filter(!is.na(entrez_gene)) %>%
    mutate(entrez_gene = as.character(entrez_gene)) %>%
    group_by(gs_name) %>%
    summarize(genes = list(unique(entrez_gene)), .groups = "drop")
  
  setNames(pathways$genes, pathways$gs_name)
}

run_fgsea_oncogenic <- function(ranks, label, pathways,
                                minSize = MIN_SIZE,
                                maxSize = MAX_SIZE,
                                nperm   = NPERM) {
  
  message("Running fgsea (MSigDB C6) for: ", label,
          " | genes=", length(ranks),
          " | sets=", length(pathways))
  
  set.seed(1)
  res <- fgsea(
    pathways = pathways,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize,
    nperm    = nperm
  ) %>%
    as_tibble() %>%
    arrange(padj) %>%
    mutate(
      # compatibilidad con tu tabla anterior
      ID          = pathway,
      Description = pathway,
      p.adjust    = padj,
      setSize     = size,
      enrichmentScore = ES,
      core_enrichment = leadingEdge %>% vapply(function(x) paste(x, collapse="/"), character(1))
    ) %>%
    select(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, core_enrichment)
  
  res
}

save_all <- function(df, label) {
  
  if (is.null(df) || nrow(df) == 0) {
    message("Sin resultados GSEA para: ", label)
    return(invisible(NULL))
  }
  
  # ---------- guardar tablas ----------
  for (od in out_dirs) {
    write_tsv(df, file.path(od, paste0("GSEA_KEGG_", label, ".tsv")))
    saveRDS(df, file.path(od, paste0("GSEA_KEGG_", label, ".rds")))
  }
  
  # ---------- plots ----------
  top_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = SHOW_TOP) %>%
    mutate(
      Direction = ifelse(NES >= 0, "Up", "Down"),
      Description = factor(Description, levels = rev(Description))
    )
  
  # 1) "Dotplot" estilo fgsea: NES vs término, tamaño=setSize, color=-log10(FDR)
  p_dot <- ggplot(top_df, aes(x = NES, y = Description)) +
    geom_point(aes(size = setSize, alpha = -log10(p.adjust))) +
    labs(
      title = paste0("GSEA (MSigDB C6 oncogenic) — ", label),
      x = "NES", y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8)
    )
  
  # 2) Barplot NES top
  p_nes <- ggplot(top_df, aes(x = Description, y = NES)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("Top oncogenic signatures — ", label),
      x = NULL, y = "NES"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8)
    )
  
  # 3) "Ridgeplot" no aplica directo sin objeto gseaResult;
  # hacemos un plot alternativo: -log10(FDR) vs término
  p_fdr <- ggplot(top_df, aes(x = -log10(p.adjust), y = Description)) +
    geom_point() +
    labs(
      title = paste0("Significance (FDR) — ", label),
      x = expression(-log[10]("FDR")), y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8)
    )
  
  # PDFs más altos para evitar solapamiento
  pdf_h <- max(10, 0.35 * nrow(top_df) + 4)  # altura dinámica
  pdf_w <- 12
  
  for (od in out_dirs) {
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_dotplot.pdf")),
           p_dot, width = pdf_w, height = pdf_h, units = "in")
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_NES_top20.pdf")),
           p_nes, width = pdf_w, height = pdf_h, units = "in")
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_ridgeplot.pdf")),
           p_fdr, width = pdf_w, height = pdf_h, units = "in")
  }
  
  invisible(df)
}

# ===================== RUN ===================== #
ranks_raw <- make_rank(DE_full_raw_DESeq2, "raw")
ranks_ae  <- make_rank(DE_full_ae_DESeq2,  "autoencoder")

pathways_onc <- get_msig_oncogenic()

df_raw <- run_fgsea_oncogenic(ranks_raw, "raw", pathways_onc)
df_ae  <- run_fgsea_oncogenic(ranks_ae,  "autoencoder", pathways_onc)

df_raw <- save_all(df_raw, "raw")
df_ae  <- save_all(df_ae,  "autoencoder")

message("DONE.")
message("Outputs overwritten in:")
message(" - ", out_dirs[1])
message(" - ", out_dirs[2])
