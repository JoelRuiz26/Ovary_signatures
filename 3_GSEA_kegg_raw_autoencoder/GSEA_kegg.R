#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(KEGGREST)
})

# ===================== CONFIG ===================== #
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

# KEGG settings
ORGANISM <- "hsa"
PVALUE_CUTOFF <- 0.05
MIN_GS <- 10
MAX_GS <- 500
EPS <- 1e-10

# Plots
SHOW_TOP <- 25
PDF_W <- 12
PDF_H_DOT <- 14
PDF_H_RIDGE <- 16
PDF_H_NES <- 14

# (Opcional) filtrar términos de infección/disease tipo COVID del REPORTE/PLOTS
FILTER_DISEASE_TERMS <- TRUE
DISEASE_PATTERN <- "COVID|CORONAVIRUS|SARS|INFLUENZA|VIRAL|INFECTION|HEPATITIS|MEASLES|TUBERCULOSIS"

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
  
  if (length(ranks) < 1000) warning("Ranking pequeño (", length(ranks), ") en ", label)
  ranks
}

# ===================== KEGG (ONE-TIME DOWNLOAD) ===================== #
# Descarga UNA vez:
# - gene -> pathway (TERM2GENE)
# - pathway -> nombre (TERM2NAME)
get_kegg_termdb_once <- function(organism = "hsa") {
  message("Descargando anotación KEGG una sola vez para: ", organism)
  
  # gene -> pathway
  # ejemplo: names = "hsa:10458", values = "path:hsa04110"
  k_link <- tryCatch({
    KEGGREST::keggLink("pathway", organism)
  }, error = function(e) {
    stop("Fallo KEGGREST::keggLink(): ", conditionMessage(e))
  })
  
  term2gene <- tibble(
    gene = sub("^hsa:", "", names(k_link)),
    term = sub("^path:", "", unname(k_link))
  ) %>%
    distinct(term, gene)
  
  # pathway -> name
  # names = "path:hsa04110", values = "Cell cycle - Homo sapiens (human)"
  k_list <- tryCatch({
    KEGGREST::keggList("pathway", organism)
  }, error = function(e) {
    stop("Fallo KEGGREST::keggList(): ", conditionMessage(e))
  })
  
  term2name <- tibble(
    term = sub("^path:", "", names(k_list)),
    name = unname(k_list)
  ) %>%
    distinct(term, name)
  
  list(term2gene = term2gene, term2name = term2name)
}

# ===================== RUN GSEA (local TERM2GENE) ===================== #
run_gsea_kegg_local <- function(ranks, label, term2gene, term2name) {
  message("Running KEGG GSEA (local TERM2GENE) for: ", label)
  
  tryCatch({
    clusterProfiler::GSEA(
      geneList     = ranks,
      TERM2GENE    = term2gene,
      TERM2NAME    = term2name,
      pvalueCutoff = PVALUE_CUTOFF,
      minGSSize    = MIN_GS,
      maxGSSize    = MAX_GS,
      eps          = EPS,
      verbose      = FALSE
    )
  }, error = function(e) {
    message("ERROR en GSEA KEGG local (", label, "): ", conditionMessage(e))
    NULL
  })
}

# ===================== FILTER (optional) ===================== #
filter_gsea_terms <- function(gsea_obj) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) return(gsea_obj)
  if (!FILTER_DISEASE_TERMS) return(gsea_obj)
  
  df <- as.data.frame(gsea_obj)
  keep <- !grepl(DISEASE_PATTERN, df$Description, ignore.case = TRUE)
  if (sum(keep) == 0) return(gsea_obj)
  
  gsea_obj@result <- gsea_obj@result[keep, , drop = FALSE]
  gsea_obj
}

# ===================== SAVE ===================== #
save_all <- function(gsea_obj, label) {
  
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("Sin resultados GSEA para: ", label)
    return(invisible(NULL))
  }
  
  gsea_obj2 <- filter_gsea_terms(gsea_obj)
  df <- as.data.frame(gsea_obj2)
  
  # ---------- guardar tablas (sobrescribe) ----------
  for (od in out_dirs) {
    write_tsv(df, file.path(od, paste0("GSEA_KEGG_", label, ".tsv")))
    saveRDS(gsea_obj2, file.path(od, paste0("GSEA_KEGG_", label, ".rds")))
  }
  
  # ---------- plots "bonitos" ----------
  p_dot <- dotplot(gsea_obj2, showCategory = SHOW_TOP, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = paste0("GSEA KEGG — ", label)) +
    theme_minimal(base_size = 12)
  
  p_ridge <- ridgeplot(gsea_obj2, showCategory = SHOW_TOP) +
    labs(title = paste0("GSEA KEGG (ridge) — ", label)) +
    theme_minimal(base_size = 12)
  
  top_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = SHOW_TOP) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  p_nes <- ggplot(top_df, aes(x = Description, y = NES)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top KEGG pathways — ", label),
         x = NULL, y = "NES") +
    theme_minimal(base_size = 12)
  
  for (od in out_dirs) {
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_dotplot.pdf")),
           p_dot, width = PDF_W, height = PDF_H_DOT, units = "in")
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_ridgeplot.pdf")),
           p_ridge, width = PDF_W, height = PDF_H_RIDGE, units = "in")
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_NES_top20.pdf")),
           p_nes, width = PDF_W, height = PDF_H_NES, units = "in")
  }
  
  invisible(df)
}

# ===================== RUN ===================== #
ranks_raw <- make_rank(DE_full_raw_DESeq2, "raw")
ranks_ae  <- make_rank(DE_full_ae_DESeq2,  "autoencoder")

# <-- KEGG se consulta UNA sola vez aquí
kegg_db <- get_kegg_termdb_once(ORGANISM)

gsea_raw <- run_gsea_kegg_local(ranks_raw, "raw", kegg_db$term2gene, kegg_db$term2name)
gsea_ae  <- run_gsea_kegg_local(ranks_ae,  "autoencoder", kegg_db$term2gene, kegg_db$term2name)

df_raw <- save_all(gsea_raw, "raw")
df_ae  <- save_all(gsea_ae,  "autoencoder")

message("DONE.")
message("Outputs overwritten in:")
message(" - ", out_dirs[1])
message(" - ", out_dirs[2])
