#!/usr/bin/env Rscript

# ============================================================
# GO over-representation of the consensus genes by direction (up / down),
# with the same universe and parameters as 4_ORA_clusterProfiler/4_ORA_core_bar.R.
#
# Input:  3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv
# Output: 9_Revision_analyses/consensus/ORA_GO_consensus_by_direction.tsv
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/consensus")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

all_genes <- read_tsv(file.path(BASE, "3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv"),
                      show_col_types = FALSE) %>%
  distinct(gene, .keep_all = TRUE)
universe  <- all_genes$gene
consensus <- all_genes %>% filter(n_sources == 3)

res <- list()
for (dir in c("up", "down")) for (ont in c("BP", "MF", "CC")) {
  g  <- consensus$gene[consensus$consensus_direction == dir]
  eg <- enrichGO(g, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = ont,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20)
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) next
  eg <- clusterProfiler::simplify(eg, cutoff = 0.3, by = "p.adjust", select_fun = min, measure = "Wang")
  res[[paste(dir, ont)]] <- as.data.frame(eg) %>%
    mutate(direction = dir, Ontology = ont, N_genes = length(g), Universe = length(universe))
}
res <- bind_rows(res) %>% arrange(direction, Ontology, p.adjust)
write_tsv(res, file.path(OUT_DIR, "ORA_GO_consensus_by_direction.tsv"))

cat("universe:", length(universe), "| up:", sum(consensus$consensus_direction == "up"),
    "down:", sum(consensus$consensus_direction == "down"), "\n")
print(res %>% group_by(direction, Ontology) %>% summarise(n = n(), .groups = "drop"))
print(res %>% group_by(direction, Ontology) %>% slice_head(n = 6) %>% ungroup() %>%
        dplyr::select(direction, Ontology, Description, GeneRatio, p.adjust), n = 60)
cat("clusterProfiler", as.character(packageVersion("clusterProfiler")),
    "org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")),
    "GO.db", as.character(packageVersion("GO.db")), "\n")
