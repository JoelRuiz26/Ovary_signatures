#!/usr/bin/env Rscript

# ============================================================
# Sensitivity of the regulator prioritization and random-gene null.
#
# Repeats the prioritization of 6_2_MR_Enriched.R / 6_3_Heatmap_TF.R
# (regulators above the |NES| elbow, enriched in highly dependent (HD)
# consensus genes, themselves dependencies and expressed) while varying
# one threshold at a time, and compares the observed result with random
# sets of HD genes drawn from the consensus.
#
# Input:  5_MRA/core_mra.tsv, 5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds,
#         5_MRA/5_0_1_Top_elbow_regulones.rds,
#         6_Depmap_ovary/6_1_0_Dependency_DEG_core.tsv,
#         6_Depmap_ovary/6_3_2_candidate_regulators.tsv (only to check the baseline),
#         7_Figures/7_2_depmap_ovary_24Q4_646_cache.rds
# Output (9_Revision_analyses/prioritization/):
#         fisher_by_regulator.tsv        Fisher test of every significant regulator
#         HD_targets_by_regulator.tsv    HD targets of each enriched regulator
#         prioritization_sensitivity.tsv one row per threshold setting
#         null_random_HD_consensus.tsv, null_random_HD_consensus_up.tsv
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

set.seed(1)
BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/prioritization")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

mra  <- read_tsv(file.path(BASE, "5_MRA/core_mra.tsv"), show_col_types = FALSE)
sigr <- mra %>% filter(padj_emp <= 0.05) %>% arrange(desc(abs(NES)))
regs <- readRDS(file.path(BASE, "5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds"))[sigr$TF]
dep  <- read_tsv(file.path(BASE, "6_Depmap_ovary/6_1_0_Dependency_DEG_core.tsv"), show_col_types = FALSE)

# Elbow of the |NES| curve, as found by 5_1_Elbow_TF.R
ELBOW <- nrow(readRDS(file.path(BASE, "5_MRA/5_0_1_Top_elbow_regulones.rds")))
cat("significant regulators:", nrow(sigr), "| elbow:", ELBOW, "\n")

# Median gene effect and expression of every gene in the ovarian lines
cache   <- readRDS(file.path(BASE, "7_Figures/7_2_depmap_ovary_24Q4_646_cache.rds"))
med_all <- cache$crispr %>% group_by(gene) %>% summarise(med = median(dep, na.rm = TRUE), .groups = "drop")
med     <- setNames(med_all$med, toupper(med_all$gene))
tpm     <- cache$expr_tpm

# Background: consensus genes with dependency data that are regulon targets.
# The dependency table has upper-case symbols (C12ORF56) and the regulons do not (C12orf56)
targets <- lapply(regs, function(x) toupper(names(x$tfmode)))
univ    <- intersect(unique(dep$gene), unique(unlist(targets)))
tg_u    <- lapply(targets, intersect, univ)
M <- lengths(tg_u)
N <- length(univ)

# Fisher test (hypergeometric tail) of every regulator against a set of HD genes
fisher_all <- function(hd) {
  hd <- intersect(hd, univ)
  n  <- length(hd)
  k  <- vapply(tg_u, function(g) sum(g %in% hd), 0L)
  p  <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  data.frame(regulon = names(tg_u), size_regulon = M, overlap = k, p_value = p, p_adj = p.adjust(p, "BH"))
}

# Prioritization: enriched regulators above the elbow, then dependency and expression filters
prioritize <- function(hd, elbow = ELBOW, min_ov = 10, alpha = 0.01, dep_thr = -0.5, tpm_thr = 1) {
  r    <- fisher_all(hd)
  top  <- sigr$TF[seq_len(min(elbow, nrow(sigr)))]
  enr  <- r$regulon[r$regulon %in% top & r$p_adj < alpha & r$overlap >= min_ov]
  cand <- enr[!is.na(med[enr]) & med[enr] <= dep_thr & !is.na(tpm[enr]) & tpm[enr] >= tpm_thr]
  list(enr = enr, cand = sort(cand))
}

hd1  <- dep$gene[dep$median_dep <= -1]
base <- prioritize(hd1)
cat("universe:", N, "| HD:", length(intersect(hd1, univ)), "| enriched:", length(base$enr),
    "| candidates:", length(base$cand), "\n", paste(base$cand, collapse = ", "), "\n")

# The baseline must reproduce the candidates of the pipeline
pipeline_cand <- fread(file.path(BASE, "6_Depmap_ovary/6_3_2_candidate_regulators.tsv"))[candidate == TRUE]$gene
if (!setequal(base$cand, pipeline_cand))
  warning("Baseline candidates differ from 6_3_2_candidate_regulators.tsv")

# Full Fisher table and HD targets of each enriched regulator
ft <- fisher_all(hd1) %>%
  left_join(sigr %>% select(regulon = TF, NES, padj_emp, rank_abs_NES), by = "regulon") %>%
  mutate(above_elbow = rank_abs_NES <= ELBOW,
         enriched = above_elbow & p_adj < 0.01 & overlap >= 10,
         median_gene_effect = med[regulon], median_tpm = tpm[regulon],
         candidate = regulon %in% base$cand) %>%
  arrange(p_adj)
write_tsv(ft, file.path(OUT_DIR, "fisher_by_regulator.tsv"))

hdt <- bind_rows(lapply(base$enr, function(r) {
  x  <- regs[[r]]
  nm <- toupper(names(x$tfmode))
  g  <- intersect(nm, intersect(hd1, univ))
  i  <- match(g, nm)
  data.frame(regulator = r, target = g, tfmode = unname(x$tfmode[i]),
             likelihood = unname(x$likelihood[i]),
             target_median_gene_effect = dep$median_dep[match(g, dep$gene)])
}))
write_tsv(hdt, file.path(OUT_DIR, "HD_targets_by_regulator.tsv"))

# Sensitivity grid: one threshold changes at a time
grid <- bind_rows(
  expand.grid(elbow = c(150, ELBOW - 30, ELBOW, ELBOW + 30, 300, 400, nrow(sigr)), min_ov = 10, alpha = 0.01, hd_thr = -1),
  expand.grid(elbow = ELBOW, min_ov = c(5, 8, 12, 15), alpha = 0.01, hd_thr = -1),
  expand.grid(elbow = ELBOW, min_ov = 10, alpha = c(0.05, 0.001), hd_thr = -1),
  expand.grid(elbow = ELBOW, min_ov = 10, alpha = 0.01, hd_thr = c(-0.5, -0.75, -1.25))) %>%
  mutate(dep_thr = -0.5, tpm_thr = 1)
# Thresholds on the regulator itself: dependency and expression
grid <- bind_rows(
  grid,
  data.frame(elbow = ELBOW, min_ov = 10, alpha = 0.01, hd_thr = -1, dep_thr = c(-0.4, -0.6, -0.75, -1), tpm_thr = 1),
  data.frame(elbow = ELBOW, min_ov = 10, alpha = 0.01, hd_thr = -1, dep_thr = -0.5, tpm_thr = c(0, 0.5, 2, 5)))
sens <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  r <- prioritize(dep$gene[dep$median_dep <= g$hd_thr], g$elbow, g$min_ov, g$alpha, g$dep_thr, g$tpm_thr)
  data.frame(g, n_HD = sum(dep$median_dep <= g$hd_thr), n_enriched = length(r$enr), n_cand = length(r$cand),
             kept_of_base = sum(base$cand %in% r$cand),
             lost = paste(setdiff(base$cand, r$cand), collapse = ","),
             gained = paste(setdiff(r$cand, base$cand), collapse = ","))
}))
print(sens)
write_tsv(sens, file.path(OUT_DIR, "prioritization_sensitivity.tsv"))

# Null: random HD genes within the consensus (all genes, and overexpressed genes only)
null_run <- function(pool, n_rep = 1000) t(replicate(n_rep, {
  r <- prioritize(sample(pool, length(intersect(hd1, univ))))
  c(enr = length(r$enr), cand = length(r$cand), base_candidates_enriched = sum(base$cand %in% r$enr))
}))
up <- intersect(dep$gene[dep$direction == "up"], univ)
for (nm in c("consensus", "consensus_up")) {
  nl <- null_run(if (nm == "consensus") univ else up)
  cat("\nNULL", nm, "(pool", if (nm == "consensus") N else length(up), "): enriched median", median(nl[, "enr"]),
      "p95", quantile(nl[, "enr"], .95), "max", max(nl[, "enr"]),
      "| candidates median", median(nl[, "cand"]), "p95", quantile(nl[, "cand"], .95), "max", max(nl[, "cand"]),
      "| P(enr >= obs) =", mean(nl[, "enr"] >= length(base$enr)),
      "| P(cand >= obs) =", mean(nl[, "cand"] >= length(base$cand)), "\n")
  write_tsv(as.data.frame(nl), file.path(OUT_DIR, paste0("null_random_HD_", nm, ".tsv")))
}

# Regulators above the elbow that pass the dependency and expression filters,
# without the enrichment filter
top   <- sigr$TF[seq_len(ELBOW)]
okdep <- top[!is.na(med[top]) & med[top] <= -0.5 & !is.na(tpm[top]) & tpm[top] >= 1]
cat("\nAbove the elbow with dependency <= -0.5 and TPM >= 1:", length(okdep),
    "| of them enriched:", sum(okdep %in% base$enr),
    "\n not enriched:", paste(setdiff(okdep, base$enr), collapse = ", "), "\n")
cat("Activated above the elbow:", sum(sigr$NES[seq_len(ELBOW)] > 0),
    "| activated among the enriched:", sum(sigr$NES[match(base$enr, sigr$TF)] > 0), "\n")
