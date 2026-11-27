#!/usr/bin/env Rscript

# ============================================================
# Fallopian tube contrast, step 4: msVIPER with the regulon of the main
# analysis and the tumor vs fallopian tube signatures (moderated t of
# each series and Stouffer z).
#
# Null: permutation of the gene labels of the signature, as
# msviper_perm_null() of 5_MRA/5_0_ovary_mra.R, with one seed per
# permutation (the result does not depend on the number of workers).
#
# Input:  5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds, 5_MRA/core_mra.tsv,
#         5_MRA/signature_core.tsv, 5_MRA/5_0_1_Top_elbow_regulones.rds,
#         6_Depmap_ovary/6_3_2_candidate_regulators.tsv,
#         9_Revision_analyses/fallopian_tube/de_all.rds, de_Stouffer.tsv (9_4_2)
# Output (9_Revision_analyses/fallopian_tube/):
#         ft_mra_results.rds, ft_mra_vs_main.tsv, ft_mra_candidates.tsv,
#         ft_mra_summary_by_contrast.tsv
# Settings: MRA_CORES parallel workers (default 10)
# ============================================================

suppressPackageStartupMessages({
  library(viper)
  library(parallel)
  library(dplyr)
})

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/fallopian_tube")
N_PERM  <- 1000; MINSIZE <- 20; SEED <- 42
N_CORES <- as.integer(Sys.getenv("MRA_CORES", "10"))

regulon <- readRDS(file.path(BASE, "5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds"))
mra  <- read.delim(file.path(BASE, "5_MRA/core_mra.tsv"))
core <- read.delim(file.path(BASE, "5_MRA/signature_core.tsv"))
cand <- read.delim(file.path(BASE, "6_Depmap_ovary/6_3_2_candidate_regulators.tsv"))
cand <- cand$gene[cand$candidate]
cat("Candidates:", length(cand), "| genes in the signature of the main analysis:", nrow(core), "\n")

# check: the regulon and the call reproduce the NES of the main analysis
o0 <- msviper(setNames(core$z_mean, core$gene), regulon, minsize = MINSIZE, verbose = FALSE)
stopifnot(setequal(names(o0$es$nes), mra$TF))
cat("Reproduction of the main analysis: max |NES difference| =", max(abs(o0$es$nes[mra$TF] - mra$NES)), "\n")

# regulators above the elbow (5_1_Elbow_TF.R): significant, ranked by |NES|
ELBOW <- nrow(readRDS(file.path(BASE, "5_MRA/5_0_1_Top_elbow_regulones.rds")))
elbow <- mra$TF[mra$padj_emp < 0.05 & mra$rank_abs_NES <= ELBOW]
stopifnot(length(elbow) == ELBOW, all(cand %in% elbow))
cat("Elbow:", ELBOW, "=", sum(mra$NES[mra$TF %in% elbow] > 0), "activated +",
    sum(mra$NES[mra$TF %in% elbow] < 0), "repressed\n")

# signatures
de <- readRDS(file.path(OUT_DIR, "de_all.rds")); st <- read.delim(file.path(OUT_DIR, "de_Stouffer.tsv"))
signatures <- lapply(de, function(d) setNames(d$t, d$gene))
signatures$Stouffer <- setNames(st$z, st$gene)
main <- c("GSE10971", "GSE69428", "Stouffer")
# sensitivity: the same signature restricted to the gene universe of the signature of the main analysis
for (nm in main) signatures[[paste0(nm, "_core_universe")]] <- signatures[[nm]][names(signatures[[nm]]) %in% core$gene]
with_perm <- c(main, paste0(main, "_core_universe"))
print(sapply(signatures, length))

permute_signature_names <- function(sig) {names(sig) <- sample(names(sig), length(sig), replace = FALSE); sig}

cl <- makeCluster(N_CORES)
clusterEvalQ(cl, {Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1"); library(viper)})
clusterExport(cl, c("regulon", "MINSIZE", "permute_signature_names"))

res <- list()
for (nm in names(signatures)) {
  sig <- signatures[[nm]]
  t0  <- Sys.time()
  obs <- msviper(sig, regulon, minsize = MINSIZE, verbose = FALSE)
  nes <- obs$es$nes; tf <- names(nes)
  out <- data.frame(TF = tf, NES = as.numeric(nes), size = as.numeric(obs$es$size[tf]),
                    p_analytic = as.numeric(obs$es$p.value[tf]), p_emp = NA_real_, padj_emp = NA_real_)
  out$padj_analytic <- p.adjust(out$p_analytic, "BH")
  if (nm %in% with_perm) {
    set.seed(SEED); perm_seeds <- sample.int(.Machine$integer.max, N_PERM)
    clusterExport(cl, c("sig", "tf", "perm_seeds"), envir = environment())
    null <- parLapply(cl, seq_len(N_PERM), function(i) {
      set.seed(perm_seeds[i])
      msviper(permute_signature_names(sig), regulon, minsize = MINSIZE, verbose = FALSE)$es$nes[tf]
    })
    null <- do.call(cbind, null)
    out$p_emp <- (1 + rowSums(abs(null) >= abs(out$NES), na.rm = TRUE)) / (N_PERM + 1)
    out$padj_emp <- p.adjust(out$p_emp, "BH")
    cat(sprintf("%s: permutation null, mean sd of the null NES = %.2f (about 1 expected if the analytic null were valid)\n",
                nm, mean(apply(null, 1, sd, na.rm = TRUE))))
  }
  res[[nm]] <- out
  cat(sprintf("%s: %d genes in the signature, %d regulators tested, %.1f min\n", nm, length(sig), nrow(out),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
saveRDS(res, file.path(OUT_DIR, "ft_mra_results.rds"))

# --- wide table against the main analysis ---
tab <- mra %>% transmute(TF, NES_main = NES, padj_main = padj_emp, rank_abs_NES,
                         above_elbow = TF %in% elbow,
                         direction_main = ifelse(NES > 0, "activated", "repressed"),
                         candidate = TF %in% cand)
for (nm in names(res)) {
  r <- res[[nm]]; i <- match(tab$TF, r$TF)
  tab[[paste0("NES_", nm)]] <- r$NES[i]
  if (nm %in% with_perm) {
    tab[[paste0("p_emp_", nm)]] <- r$p_emp[i]; tab[[paste0("padj_emp_", nm)]] <- r$padj_emp[i]
  }
  tab[[paste0("padj_analytic_", nm)]] <- r$padj_analytic[i]
}
tab <- tab[order(tab$rank_abs_NES), ]
write.table(tab, file.path(OUT_DIR, "ft_mra_vs_main.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- candidates ---
ct <- tab[tab$candidate, c("TF", "NES_main", "padj_main", "rank_abs_NES",
                           unlist(lapply(main, function(n) paste0(c("NES_", "p_emp_", "padj_emp_"), n))))]
write.table(ct, file.path(OUT_DIR, "ft_mra_candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250); print(ct %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)

# --- summary by contrast ---
smry <- list()
for (nm in names(res)) {
  nes  <- tab[[paste0("NES_", nm)]]
  padj <- if (nm %in% with_perm) tab[[paste0("padj_emp_", nm)]] else tab[[paste0("padj_analytic_", nm)]]
  ok   <- !is.na(nes)
  sig_same <- sign(nes) == sign(tab$NES_main) & padj < 0.05
  for (set in c("all_common", "main_significant", "above_elbow", "elbow_activated", "elbow_repressed", "candidates")) {
    k <- ok & switch(set, all_common = TRUE, main_significant = tab$padj_main < 0.05,
                     above_elbow = tab$above_elbow,
                     elbow_activated = tab$above_elbow & tab$NES_main > 0,
                     elbow_repressed = tab$above_elbow & tab$NES_main < 0, candidates = tab$candidate)
    ms <- sum(sign(nes[k]) == sign(tab$NES_main[k]))
    smry[[length(smry) + 1]] <- data.frame(
      contrast = nm, p_type = ifelse(nm %in% with_perm, "gene_permutation", "msviper_analytic"), set = set,
      n_in_main = switch(set, all_common = nrow(tab), main_significant = sum(tab$padj_main < 0.05),
                         above_elbow = ELBOW, elbow_activated = sum(tab$above_elbow & tab$NES_main > 0),
                         elbow_repressed = sum(tab$above_elbow & tab$NES_main < 0), candidates = length(cand)),
      n_evaluated = sum(k),
      rho_spearman = cor(nes[k], tab$NES_main[k], method = "spearman"),
      n_same_sign = ms, frac_same_sign = ms / sum(k),
      p_binomial = binom.test(ms, sum(k), 0.5)$p.value,
      n_sig_same_sign = sum(sig_same[k]), frac_sig_same_sign = mean(sig_same[k]),
      n_sig_opposite_sign = sum(sign(nes[k]) != sign(tab$NES_main[k]) & padj[k] < 0.05))
  }
}
smry <- do.call(rbind, smry)
write.table(smry, file.path(OUT_DIR, "ft_mra_summary_by_contrast.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(smry %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)

# selected repressed regulators (also among the markers printed by 9_4_2)
cat("\nSelected repressed regulators:\n")
print(tab[tab$TF %in% c("ZEB1", "PRDM6", "OSR1", "LHX9"),
          c("TF", "NES_main", "rank_abs_NES", unlist(lapply(main, function(n) paste0(c("NES_", "padj_emp_"), n))))] %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
sessionInfo()
