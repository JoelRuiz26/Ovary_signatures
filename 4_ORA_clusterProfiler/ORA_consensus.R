#!/usr/bin/env Rscript

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(enrichplot)
  library(ggplot2)
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

deseq_sig_signed <- function(df, source_name) {
  df %>%
    mutate(
      gene = clean_gene(Symbol),
      direction = case_when(
        !is.na(pvalue) & pvalue < P_CUTOFF &
          !is.na(log2FoldChange) & log2FoldChange >  LFC_CUTOFF ~ "up",
        !is.na(pvalue) & pvalue < P_CUTOFF &
          !is.na(log2FoldChange) & log2FoldChange < -LFC_CUTOFF ~ "down",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(gene), !is.na(direction)) %>%
    distinct(gene, direction) %>%
    mutate(source = source_name)
}

geo_sig_signed <- function(df) {
  df %>%
    mutate(
      gene = clean_gene(gene),
      direction = tolower(trimws(direction))
    ) %>%
    filter(
      !is.na(gene),
      !is.na(p_raw), p_raw < P_CUTOFF,
      direction %in% c("up", "down")
    ) %>%
    distinct(gene, direction) %>%
    mutate(source = "GEO")
}

# ============================================================
# LOAD DATA
# ============================================================
DE_raw <- read_rds_safe(paths$raw_deseq2)
DE_ae  <- read_rds_safe(paths$ae_deseq2)
GEO    <- read_rds_safe(paths$geo_rra)

raw_sig <- deseq_sig_signed(DE_raw, "TCGA_GTEx")
ae_sig  <- deseq_sig_signed(DE_ae,  "TCGA_AE")
geo_sig <- geo_sig_signed(GEO)

# ============================================================
# SPLIT BY DIRECTION
# ============================================================
raw_up   <- raw_sig %>% filter(direction == "up")   %>% pull(gene)
raw_down <- raw_sig %>% filter(direction == "down") %>% pull(gene)

ae_up    <- ae_sig  %>% filter(direction == "up")   %>% pull(gene)
ae_down  <- ae_sig  %>% filter(direction == "down") %>% pull(gene)

geo_up   <- geo_sig %>% filter(direction == "up")   %>% pull(gene)
geo_down <- geo_sig %>% filter(direction == "down") %>% pull(gene)

# ============================================================
# CONSENSUS (3/3) + UNIVERSE (UNION) PER DIRECTION
# ============================================================
up_consensus   <- Reduce(intersect, list(raw_up, ae_up, geo_up))
down_consensus <- Reduce(intersect, list(raw_down, ae_down, geo_down))

up_universe    <- Reduce(union, list(raw_up, ae_up, geo_up))
down_universe  <- Reduce(union, list(raw_down, ae_down, geo_down))

cat("[INFO] UP consensus 3/3:", length(up_consensus), "\n")
cat("[INFO] DOWN consensus 3/3:", length(down_consensus), "\n")
cat("[INFO] UP universe:", length(up_universe), "\n")
cat("[INFO] DOWN universe:", length(down_universe), "\n")

# ============================================================
# ORA FUNCTION (MISMA LÓGICA QUE TU PIPELINE)
# ============================================================
perform_ora <- function(genes, universe, ontology, direction_label) {
  
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
    mutate(Direction = direction_label,
           Ontology  = ontology,
           N_genes   = length(genes),
           Universe  = length(universe))
}

# ============================================================
# RUN ORA (UP & DOWN × BP/MF/CC)
# ============================================================
ora_results <- bind_rows(
  lapply(ONTOLOGIES, function(ont) {
    bind_rows(
      perform_ora(up_consensus,   up_universe,   ont, "UP"),
      perform_ora(down_consensus, down_universe, ont, "DOWN")
    )
  })
)

# ============================================================
# SAVE RESULTS
# ============================================================
out_tsv <- file.path(
  out_dir,
  paste0("ORA_GO_consensus3of3_UP_DOWN_BP_MF_CC_SIMPLIFIED_", SIMPLIFY_CUTOFF, ".tsv")
)

write_tsv(ora_results, out_tsv)

cat("ORA results saved to:\n", out_tsv, "\n")

# ============================================================
# QUICK SUMMARY
# ============================================================
summary_df <- ora_results %>%
  count(Direction, Ontology, name = "n_terms")

print(summary_df)

end_time <- Sys.time()
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

save.image("~/Ovary_signatures/4_ORA_clusterProfiler/ORA_IMAGE.RData")

