#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
  library(stringr)
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
MSIG_COLLECTION <- "C6"   # oncogenic signatures
MIN_SIZE <- 15
MAX_SIZE <- 500
SHOW_TOP <- 30

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

# ---- MSigDB C6: en tu versión, Entrez está en ncbi_gene ----
get_msig_oncogenic <- function() {
  message("Cargando MSigDB ", MSIG_COLLECTION, " (oncogenic signatures)...")
  
  msig <- msigdbr(species = "Homo sapiens", collection = MSIG_COLLECTION)
  
  if (!("ncbi_gene" %in% colnames(msig))) {
    stop("No encontré 'ncbi_gene' en msigdbr. Columnas:\n", paste(colnames(msig), collapse = ", "))
  }
  
  # mapping para labels legibles (puede repetirse, por eso luego haremos Label único)
  map_tbl <- msig %>%
    distinct(gs_name, gs_description) %>%
    mutate(gs_description = ifelse(is.na(gs_description) | gs_description == "", gs_name, gs_description))
  
  pathways_tbl <- msig %>%
    transmute(
      gs_name = gs_name,
      entrez  = as.character(ncbi_gene)
    ) %>%
    filter(!is.na(entrez), entrez != "") %>%
    group_by(gs_name) %>%
    summarize(genes = list(unique(entrez)), .groups = "drop")
  
  pathways <- setNames(pathways_tbl$genes, pathways_tbl$gs_name)
  list(pathways = pathways, map_tbl = map_tbl)
}

run_fgsea_oncogenic <- function(ranks, label, pathways,
                                minSize = MIN_SIZE,
                                maxSize = MAX_SIZE) {
  message("Running fgseaMultilevel (MSigDB C6) for: ", label,
          " | genes=", length(ranks),
          " | sets=", length(pathways))
  
  set.seed(1)
  fgseaMultilevel(
    pathways = pathways,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize
  ) %>%
    as_tibble() %>%
    arrange(padj) %>%
    mutate(
      ID = pathway,
      p.adjust = padj,
      pvalue   = pval,
      setSize  = size,
      enrichmentScore = ES,
      core_enrichment = vapply(leadingEdge, function(x) paste(x, collapse="/"), character(1))
    ) %>%
    select(ID, setSize, enrichmentScore, NES, pvalue, p.adjust, core_enrichment)
}

save_all <- function(df, label, map_tbl) {
  
  if (is.null(df) || nrow(df) == 0) {
    message("Sin resultados GSEA para: ", label)
    return(invisible(NULL))
  }
  
  # añade descripción bonita
  df2 <- df %>%
    left_join(map_tbl, by = c("ID" = "gs_name")) %>%
    mutate(
      Description = ifelse(is.na(gs_description) | gs_description == "", ID, gs_description),
      Description = str_wrap(Description, width = 55)
    ) %>%
    select(ID, Description, everything(), -gs_description)
  
  # ---------- guardar tablas (sobrescribe) ----------
  for (od in out_dirs) {
    write_tsv(df2, file.path(od, paste0("GSEA_KEGG_", label, ".tsv")))
    saveRDS(df2, file.path(od, paste0("GSEA_KEGG_", label, ".rds")))
  }
  
  # ---------- plots ----------
  # FIX: evitar factor levels duplicados creando un Label único por fila
  top_df <- df2 %>%
    arrange(p.adjust) %>%
    slice_head(n = SHOW_TOP) %>%
    mutate(
      Label = paste0(Description, " (", ID, ")"),
      Label = factor(Label, levels = rev(unique(Label)))
    )
  
  pdf_h <- max(14, 0.42 * nrow(top_df) + 6)
  pdf_w <- 13
  
  # 1) Dotplot: color por NES, size por setSize, alpha por FDR
  p_dot <- ggplot(top_df, aes(x = NES, y = Label)) +
    geom_point(aes(size = setSize, color = NES, alpha = -log10(p.adjust))) +
    scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
    guides(alpha = "none") +
    labs(
      title = paste0("GSEA (MSigDB C6 oncogenic) — ", label),
      x = "NES", y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))
  
  # 2) NES barplot con signo
  p_nes <- ggplot(top_df, aes(x = Label, y = NES, fill = NES)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
    guides(fill = "none") +
    labs(
      title = paste0("Top oncogenic signatures — ", label),
      x = NULL, y = "NES"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))
  
  # 3) Significance: -log10(FDR)
  p_fdr <- ggplot(top_df, aes(x = -log10(p.adjust), y = Label)) +
    geom_point(aes(size = setSize, color = NES), alpha = 0.9) +
    scale_color_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
    labs(
      title = paste0("Significance (FDR) — ", label),
      x = expression(-log[10]("FDR")), y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))
  
  for (od in out_dirs) {
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_dotplot.pdf")),
           p_dot, width = pdf_w, height = pdf_h, units = "in")
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_NES_top20.pdf")),
           p_nes, width = pdf_w, height = pdf_h, units = "in")
    ggsave(file.path(od, paste0("PLOT_GSEA_KEGG_", label, "_ridgeplot.pdf")),
           p_fdr, width = pdf_w, height = pdf_h, units = "in")
  }
  
  invisible(df2)
}

# ===================== RUN ===================== #
ranks_raw <- make_rank(DE_full_raw_DESeq2, "raw")
ranks_ae  <- make_rank(DE_full_ae_DESeq2,  "autoencoder")

msig_obj <- get_msig_oncogenic()
pathways_onc <- msig_obj$pathways
map_tbl <- msig_obj$map_tbl

df_raw <- run_fgsea_oncogenic(ranks_raw, "raw", pathways_onc)
df_ae  <- run_fgsea_oncogenic(ranks_ae,  "autoencoder", pathways_onc)

df_raw <- save_all(df_raw, "raw", map_tbl)
df_ae  <- save_all(df_ae,  "autoencoder", map_tbl)

message("DONE.")
message("Outputs overwritten in:")
message(" - ", out_dirs[1])
message(" - ", out_dirs[2])
