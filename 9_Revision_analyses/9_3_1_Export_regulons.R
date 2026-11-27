#!/usr/bin/env Rscript

# ============================================================
# Regulons (target, mode of regulation, likelihood) of the enriched
# regulators evaluated in 6_3_Heatmap_TF.R, flagged as candidate or not,
# as a flat table for the Perturb-seq analysis (9_3_2).
#
# Input:  5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds,
#         6_Depmap_ovary/6_3_2_candidate_regulators.tsv
# Output: 9_Revision_analyses/perturbseq/enriched_regulons.tsv
# ============================================================

suppressPackageStartupMessages(library(data.table))

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/perturbseq")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

reg  <- readRDS(file.path(BASE, "5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds"))
cand <- fread(file.path(BASE, "6_Depmap_ovary/6_3_2_candidate_regulators.tsv"))

r <- rbindlist(lapply(cand$gene, function(g) data.table(
  driver = g, candidate = cand[gene == g]$candidate, target = names(reg[[g]]$tfmode),
  tfmode = reg[[g]]$tfmode, likelihood = reg[[g]]$likelihood)))
fwrite(r, file.path(OUT_DIR, "enriched_regulons.tsv"), sep = "\t")
print(r[candidate == TRUE, .(n = .N, pos = sum(tfmode > 0)), by = driver])
