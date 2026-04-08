#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
  library(scales)
})

options(stringsAsFactors = FALSE)

# ============================================================
# CONFIG
# ============================================================
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw_deseq2 = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae_deseq2  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  geo_rra    = file.path(base_dir, "2_DEG_GEO",  "GEO_ovarian_cancer_RRA_results.rds")
)

out_dir <- file.path(base_dir, "4_ORA_clusterProfiler")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# thresholds (significance definition)
P_CUTOFF   <- 0.01
LFC_CUTOFF <- 1

# ORA thresholds
P_ORA <- 0.05
Q_ORA <- 0.20

# semantic reduction (MISMO CRITERIO QUE TU SCRIPT)
SIMPLIFY_CUTOFF  <- 0.3
SIMPLIFY_BY      <- "p.adjust"
SIMPLIFY_MEASURE <- "Wang"
SIMPLIFY_SELECT  <- min

ONTOLOGIES <- c("BP", "MF", "CC")

# ============================================================
# HELPERS
# ============================================================
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("File not found: ", p)
  readRDS(p)
}

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# ============================================================
# SIGNIFICANT GENES: ALL TOGETHER (NO UP/DOWN SPLIT)
# ============================================================
deseq_sig_all <- function(df, source_name) {
  df %>%
    mutate(
      gene = clean_gene(Symbol),
      is_sig = !is.na(pvalue) & pvalue < P_CUTOFF &
        !is.na(log2FoldChange) & abs(log2FoldChange) > LFC_CUTOFF
    ) %>%
    filter(!is.na(gene), is_sig) %>%
    distinct(gene) %>%
    mutate(source = source_name)
}

geo_sig_all <- function(df) {
  df %>%
    mutate(
      gene = clean_gene(gene)
    ) %>%
    filter(
      !is.na(gene),
      !is.na(p_raw), p_raw < P_CUTOFF
    ) %>%
    distinct(gene) %>%
    mutate(source = "GEO")
}

# ============================================================
# LOAD DATA
# ============================================================
DE_raw <- read_rds_safe(paths$raw_deseq2)
DE_ae  <- read_rds_safe(paths$ae_deseq2)
GEO    <- read_rds_safe(paths$geo_rra)

raw_sig <- deseq_sig_all(DE_raw, "TCGA_GTEx")
ae_sig  <- deseq_sig_all(DE_ae,  "TCGA_AE")
geo_sig <- geo_sig_all(GEO)

# ============================================================
# GENE LISTS (ALL SIGNIFICANT GENES)
# ============================================================
raw_genes <- raw_sig %>% pull(gene)
ae_genes  <- ae_sig  %>% pull(gene)
geo_genes <- geo_sig %>% pull(gene)

# ============================================================
# CONSENSUS (3/3) + UNIVERSE (UNION)
# ============================================================
consensus_genes <- Reduce(intersect, list(raw_genes, ae_genes, geo_genes))
universe_genes   <- Reduce(union, list(raw_genes, ae_genes, geo_genes))

cat("[INFO] Consensus 3/3:", length(consensus_genes), "\n")
cat("[INFO] Universe:", length(universe_genes), "\n")

# ============================================================
# ORA FUNCTION
# ============================================================
perform_ora <- function(genes, universe, ontology) {
  
  if (length(genes) == 0) return(data.frame())
  
  eg <- tryCatch({
    enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = P_ORA,
      qvalueCutoff  = Q_ORA,
      readable      = FALSE
    )
  }, error = function(e) NULL)
  
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) return(data.frame())
  
  # pairwise similarity (necesario para simplify robusto)
  eg <- tryCatch({
    enrichplot::pairwise_termsim(eg, method = SIMPLIFY_MEASURE)
  }, error = function(e) eg)
  
  eg <- tryCatch({
    clusterProfiler::simplify(
      eg,
      cutoff     = SIMPLIFY_CUTOFF,
      by         = SIMPLIFY_BY,
      select_fun = SIMPLIFY_SELECT,
      measure    = SIMPLIFY_MEASURE
    )
  }, error = function(e) eg)
  
  df <- as.data.frame(eg)
  if (nrow(df) == 0) return(df)
  
  df %>%
    arrange(p.adjust, desc(Count)) %>%
    mutate(
      Group    = "All genes",
      Ontology = ontology,
      N_genes  = length(genes),
      Universe = length(universe)
    )
}

# ============================================================
# RUN ORA (ALL TOGETHER × BP/MF/CC)
# ============================================================
ora_results <- bind_rows(
  lapply(ONTOLOGIES, function(ont) {
    perform_ora(consensus_genes, universe_genes, ont)
  })
)

# ============================================================
# SAVE RESULTS
# ============================================================
out_tsv <- file.path(
  out_dir,
  paste0("ORA_GO_consensus3of3_ALL_BP_MF_CC_SIMPLIFIED_", SIMPLIFY_CUTOFF, ".tsv")
)

write_tsv(ora_results, out_tsv)

cat("ORA results saved to:\n", out_tsv, "\n")

# ============================================================
# DOTPLOT COMPACTO Y ELEGANTE
# ============================================================
parse_ratio <- function(x) {
  sapply(x, function(z) {
    if (is.na(z) || !nzchar(z)) return(NA_real_)
    parts <- strsplit(z, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_real_)
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

TOP_N_PER_ONTOLOGY <- 12

plot_df <- ora_results %>%
  mutate(
    Ontology  = factor(Ontology, levels = c("BP", "MF", "CC")),
    GeneRatio_num = parse_ratio(GeneRatio),
    mlog10FDR = -log10(p.adjust),
    Description_wrapped = str_wrap(Description, width = 36)
  ) %>%
  filter(!is.na(p.adjust), !is.na(Count)) %>%
  group_by(Ontology) %>%
  slice_min(order_by = p.adjust, n = TOP_N_PER_ONTOLOGY, with_ties = FALSE) %>%
  ungroup()

# ordenar términos dentro de cada panel
plot_df <- plot_df %>%
  group_by(Ontology) %>%
  mutate(
    Description_wrapped = factor(
      Description_wrapped,
      levels = rev(unique(Description_wrapped[order(p.adjust)]))
    )
  ) %>%
  ungroup()

p_dot <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description_wrapped)) +
  geom_col(
    aes(fill = mlog10FDR),
    width = 0.75
  ) +
  scale_fill_viridis_c(
    option = "D",
    direction = 1,
    limits = c(2, 10),
    breaks = c(2, 4, 6, 8, 10),
    labels = c("2", "4", "6", "8", "10"),
    oob = scales::squish,
    name = expression(-log[10]~"(FDR)")
  ) +
  scale_size(
    range = c(2.5, 8),
    name = "Gene count"
  ) +
  facet_grid(
    Ontology ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  labs(
    title = "GO ORA (top terms; semantic-reduced)",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14) +   # 🔥 AUMENTA TODO GLOBAL
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 18),
    
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(face = "bold", size = 12),
    
    axis.title.x = element_text(face = "bold", size = 14),
    
    strip.text.y = element_text(face = "bold", size = 13, angle = 0),
    
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12),
    
    legend.key.height = unit(0.7, "cm")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(4.5, "cm"),
      barwidth  = unit(0.6, "cm"),
      ticks.colour = "grey20",
      frame.colour = "grey40"
    )
  )

# ============================================================
# SAVE PLOT
# ============================================================
ggsave(
  file.path(out_dir, "ORA_GO_topTerms_dotplot_COMPACT_ALL.png"),
  p_dot,
  width  = 7.8,
  height = 8.8,
  dpi    = 600
)

ggsave(
  file.path(out_dir, "ORA_GO_topTerms_dotplot_COMPACT_ALL.pdf"),
  p_dot,
  width  = 7.8,
  height = 8.8
)

cat("Compact dotplot saved.\n")

# ============================================================
# QUICK SUMMARY
# ============================================================
summary_df <- ora_results %>%
  group_by(Ontology) %>%
  summarise(n_terms = n(), .groups = "drop")

print(summary_df)

end_time <- Sys.time()
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

