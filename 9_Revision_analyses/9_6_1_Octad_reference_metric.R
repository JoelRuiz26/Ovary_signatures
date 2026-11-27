#!/usr/bin/env Rscript

# ============================================================
# Similarity metric used by octad::computeRefTissue to select the
# reference samples of the autoencoder (AE) contrast.
#
# Runs the call of 1_DGE_AE/1_1_Deseq_octad_AE.R, recomputes the metric
# step by step (median Spearman correlation of each candidate normal or
# adjacent sample with the TCGA-OV tumors, on the precomputed autoencoder
# features of OCTAD) and annotates every candidate sample.
#
# Input:  octad / octad.db (ExperimentHub EH7274 and EH7265);
#         3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_1_expr_raw_Tumor_Ctl*.rds
#         (optional, to check the samples against those used by the pipeline)
# Output (9_Revision_analyses/octad_metric/):
#         metric_by_sample.tsv, metric_by_tissue.tsv, group_summary.tsv,
#         composition_AE_controls.tsv, composition_AE_controls_by_group.tsv,
#         tumor_vs_tumor_context.tsv,
#         octad_out/case_normal_median_cor.csv (written by the package)
# ============================================================

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
})

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
out     <- file.path(BASE, "9_Revision_analyses/octad_metric")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
CONTROL_SIZE <- 88   # as in 1_1_Deseq_octad_AE.R

phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
expSet  <- octad.db::get_ExperimentHub_data("EH7265")

# Cases and GTEx controls: same filters as the pipeline scripts
case_ovary <- phenoDF %>%
  filter(grepl("OVARY", biopsy.site, ignore.case = TRUE), data.source == "TCGA", sample.type != "adjacent") %>%
  pull(sample.id)
ctrl_gtex <- phenoDF %>%
  filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
         (data.source == "GTEX" & sample.type == "normal") | (data.source == "TCGA" & sample.type == "adjacent")) %>%
  pull(sample.id)
cat("Cases:", length(case_ovary), "| GTEx ovary controls:", length(ctrl_gtex), "\n")
print(table(phenoDF$data.source[match(ctrl_gtex, phenoDF$sample.id)], phenoDF$sample.type[match(ctrl_gtex, phenoDF$sample.id)]))

# 1) The call of the pipeline, twice with different seeds (determinism check)
run_ref <- function(seed, folder) {
  set.seed(seed)
  dir.create(folder, showWarnings = FALSE)
  suppressWarnings(suppressMessages(
    computeRefTissue(case_id = case_ovary, adjacent = TRUE, source = "octad", control_size = CONTROL_SIZE, outputFolder = folder)
  )) %>% as.character()
}
sel1 <- run_ref(1, file.path(out, "octad_out"))
sel2 <- run_ref(20261126, file.path(out, "octad_out_seed2"))
cat("Selected:", length(sel1), "| identical and in the same order with another seed:", identical(sel1, sel2), "\n")

# 2) Independent recomputation, step by step
normal_id <- c(phenoDF$sample.id[phenoDF$sample.type == "adjacent"], phenoDF$sample.id[phenoDF$sample.type == "normal"])
normal_id <- normal_id[normal_id %in% colnames(expSet)]
cases <- case_ovary[case_ovary %in% colnames(expSet)]
cat("Candidates:", length(normal_id), "| cases in the matrix:", length(cases), "| features:", nrow(expSet), "\n")
iqr <- apply(expSet[, normal_id], 1, stats::IQR)
vg  <- order(iqr, decreasing = TRUE)[seq_len(min(500, length(iqr)))]
cat("Features kept by the IQR filter:", length(vg), "of", nrow(expSet), "\n")
rho <- cor(expSet[vg, normal_id], expSet[vg, cases], method = "spearman")
med <- apply(rho, 1, median)
# the result is the same without the IQR step (all features pass; Spearman does not depend on the row order)
med_noiqr <- apply(cor(expSet[, normal_id], expSet[, cases], method = "spearman"), 1, median)
cat("Max difference with/without the IQR step:", max(abs(med - med_noiqr)), "\n")

# comparison with the file written by the package
csv <- read.csv(file.path(out, "octad_out", "case_normal_median_cor.csv"))
cat("Max difference vs case_normal_median_cor.csv of the package:", max(abs(csv$cor - med[csv$sample.id])),
    "| same samples:", setequal(csv$sample.id, names(med)), "\n")
cat("The selected samples are the first", CONTROL_SIZE, "of the package csv:", identical(sel1, csv$sample.id[seq_len(CONTROL_SIZE)]), "\n")
cat("Quantile of the cut-off (cor_cutoff = '0' -> 0%):", min(med), "-> no candidate is dropped by the threshold\n")

# 3) Comparison with the samples used by the pipeline (3_0_Counts_norm.R), when available
f_ae   <- file.path(BASE, "3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_1_expr_raw_Tumor_CtlAuto.rds")
f_gtex <- file.path(BASE, "3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_1_expr_raw_Tumor_Ctlhomolog.rds")
if (file.exists(f_ae) && file.exists(f_gtex)) {
  ids_ae <- colnames(readRDS(f_ae)); ids_gtex <- colnames(readRDS(f_gtex))
  cat("AE controls = controls of the pipeline matrix:", setequal(sel1, setdiff(ids_ae, case_ovary)), "\n")
  cat("GTEx controls = controls of the pipeline matrix:", setequal(ctrl_gtex, setdiff(ids_gtex, case_ovary)), "\n")
  cat("Cases = cases of the pipeline matrices:", all(case_ovary %in% ids_ae), all(case_ovary %in% ids_gtex), "\n")
} else cat("Pipeline matrices of 3_0_Counts_norm.R not found: comparison skipped\n")

# 4) Table by sample; rank 1 = highest correlation, in the order used by the function
ord <- csv$sample.id
group_tissue <- function(x) {
  x <- toupper(x)
  x <- sub("^BRAIN.*", "BRAIN", x); x <- sub("^KIDNEY.*", "KIDNEY", x); x <- sub("^THYROID.*", "THYROID", x)
  x <- sub("^BREAST.*", "BREAST", x); x <- sub("^(UTERUS|ENDOMETRIUM)$", "UTERUS/ENDOMETRIUM", x)
  x
}
tab <- tibble(sample.id = ord, median_spearman = unname(med[ord]), rank = seq_along(ord)) %>%
  left_join(phenoDF %>% select(sample.id, data.source, sample.type, biopsy.site, gender), by = "sample.id") %>%
  mutate(rank_min = rank(-median_spearman, ties.method = "min"),
         percentile = 100 * (1 - (rank - 1) / n()),
         source = case_when(data.source == "GTEX" ~ "GTEx", data.source == "TCGA" ~ "TCGA adjacent", data.source == "TARGET" ~ "TARGET adjacent", TRUE ~ data.source),
         tissue = biopsy.site, tissue_group = group_tissue(biopsy.site),
         selected_AE = sample.id %in% sel1, ovary_GTEx = sample.id %in% ctrl_gtex) %>%
  select(sample.id, median_spearman, rank, rank_min, percentile, source, sample.type, tissue, tissue_group, gender, selected_AE, ovary_GTEx)
stopifnot(!anyNA(tab$source), sum(tab$selected_AE) == CONTROL_SIZE, sum(tab$ovary_GTEx) == length(ctrl_gtex))
write.table(tab, file.path(out, "metric_by_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 5) Key figures
cat("\nTies: repeated values of the metric:", sum(duplicated(tab$median_spearman)),
    "| value at the cut-off rank:", tab$median_spearman[CONTROL_SIZE], "| next value:", tab$median_spearman[CONTROL_SIZE + 1], "\n")
smry <- function(x) c(n = length(x), min = min(x), q25 = unname(quantile(x, .25)), median = median(x), q75 = unname(quantile(x, .75)), max = max(x))
groups <- list(AE_selected = tab$median_spearman[tab$selected_AE], ovary_GTEx = tab$median_spearman[tab$ovary_GTEx],
               rest_not_AE = tab$median_spearman[!tab$selected_AE],
               rest_not_AE_not_ovary = tab$median_spearman[!tab$selected_AE & !tab$ovary_GTEx],
               all = tab$median_spearman)
num <- do.call(rbind, lapply(groups, smry)); print(round(num, 4))
write.table(data.frame(group = rownames(num), num), file.path(out, "group_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
shared <- intersect(sel1, ctrl_gtex); cat("\nShared, AE and GTEx ovary controls:", length(shared), ":", shared, "\n")
print(tab %>% filter(sample.id %in% shared) %>% select(sample.id, median_spearman, rank))
cat("\nRanks of the GTEx ovary controls:\n"); print(summary(tab$rank[tab$ovary_GTEx]))
tops <- c(CONTROL_SIZE, 100, 200, 500, 1000)
cat("GTEx ovary within the top", paste(tops, collapse = "/"), ":", sapply(tops, function(k) sum(tab$rank[tab$ovary_GTEx] <= k)), "\n")
cat("\nComposition of the selected samples (tissue x source):\n")
c1 <- tab %>% filter(selected_AE) %>% dplyr::count(tissue_group, tissue, source, name = "n") %>% arrange(desc(n)); print(as.data.frame(c1))
write.table(c1, file.path(out, "composition_AE_controls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
c2 <- tab %>% filter(selected_AE) %>% group_by(tissue_group) %>%
  summarise(n = n(), min = min(median_spearman), median = median(median_spearman), max = max(median_spearman), best_rank = min(rank)) %>% arrange(desc(n))
print(as.data.frame(c2)); write.table(c2, file.path(out, "composition_AE_controls_by_group.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSelected samples by source:\n"); print(table(tab$source[tab$selected_AE]))
cat("\nCandidates by source:\n"); print(table(tab$source))
# summary by tissue of all the candidates
pt <- tab %>% group_by(tissue_group) %>%
  summarise(n = n(), n_AE_selected = sum(selected_AE), min = min(median_spearman), median = median(median_spearman), max = max(median_spearman), best_rank = min(rank)) %>% arrange(desc(median))
write.table(pt, file.path(out, "metric_by_tissue.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nTissues (all candidates), ordered by median:\n"); print(as.data.frame(pt), digits = 3)

# 6) Context: the same metric among tumors (each tumor against the others)
tt <- cor(expSet[, cases], method = "spearman"); diag(tt) <- NA
mt <- apply(tt, 1, median, na.rm = TRUE)
cat("\nTumor vs tumors (median Spearman per tumor): "); print(round(smry(mt), 4))
write.table(data.frame(sample.id = names(mt), median_spearman_vs_other_tumors = mt), file.path(out, "tumor_vs_tumor_context.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# spread within each row: IQR of the correlation of each selected sample with the tumors
cat("IQR per selected sample of its correlation with the tumors: "); print(round(summary(apply(rho[sel1, ], 1, IQR)), 4))
sessionInfo()
