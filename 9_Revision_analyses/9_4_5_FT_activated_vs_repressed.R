#!/usr/bin/env Rscript

# ============================================================
# Fallopian tube contrast, step 5: activated against repressed
# regulators when the reference is fallopian tube epithelium, and
# context figures for the text.
#
# Input:  9_Revision_analyses/fallopian_tube/ft_mra_vs_main.tsv (9_4_4),
#         consensus_concordance.tsv (9_4_3), de_all.rds, de_Stouffer.tsv (9_4_2),
#         5_MRA/signature_core.tsv
# Output (9_Revision_analyses/fallopian_tube/):
#         elbow_activated_vs_repressed.tsv,
#         ft_mra_vs_main_ranked.tsv (ft_mra_vs_main.tsv plus the rank by |NES|
#         of each regulator within each fallopian tube analysis)
# ============================================================

suppressPackageStartupMessages(library(dplyr))
options(width = 220)

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/fallopian_tube")

tab  <- read.delim(file.path(OUT_DIR, "ft_mra_vs_main.tsv"))
con  <- read.delim(file.path(OUT_DIR, "consensus_concordance.tsv"))
main <- c("GSE10971", "GSE69428", "Stouffer")
ELBOW <- sum(tab$above_elbow)

cat("== Consensus: up vs down (Fisher's exact test on the same-sign counts) ==\n")
for (nm in main) {
  u <- con[con$contrast == nm & con$consensus_direction == "up", ]
  d <- con[con$contrast == nm & con$consensus_direction == "down", ]
  ft <- fisher.test(matrix(c(u$n_same_sign, u$n_measured - u$n_same_sign, d$n_same_sign, d$n_measured - d$n_same_sign), 2))
  cat(sprintf("%s: up %d/%d (%.1f%%) vs down %d/%d (%.1f%%); OR = %.2f, p = %.2g\n", nm, u$n_same_sign, u$n_measured,
              100 * u$frac_same_sign, d$n_same_sign, d$n_measured, 100 * d$frac_same_sign, ft$estimate, ft$p.value))
}

cat("\n== Baseline: regulators significant with the fallopian tube signature (all those shared with the main analysis) ==\n")
rows <- list()
for (nm in main) {
  nes <- tab[[paste0("NES_", nm)]]; pe <- tab[[paste0("padj_emp_", nm)]]; pa <- tab[[paste0("padj_analytic_", nm)]]
  ok <- !is.na(nes); elbow <- tab$above_elbow & ok
  same <- sign(nes) == sign(tab$NES_main)
  cat(sprintf("%s: %d shared; significant by permutation %d (%.1f%%), analytic %d (%.1f%%)\n", nm, sum(ok),
              sum(pe[ok] < 0.05), 100 * mean(pe[ok] < 0.05), sum(pa[ok] < 0.05), 100 * mean(pa[ok] < 0.05)))
  # rank by |NES| within the fallopian tube analysis (shared regulators)
  rk <- rep(NA_integer_, nrow(tab)); rk[ok] <- rank(-abs(nes[ok]), ties.method = "first")
  tab[[paste0("rank_abs_NES_", nm)]] <- rk
  for (s in c("activated", "repressed")) {
    k <- elbow & tab$direction_main == s
    rows[[length(rows) + 1]] <- data.frame(
      contrast = nm, direction = s, n = sum(k),
      same_sign = sum(same[k]),
      sig_perm_same_sign = sum(same[k] & pe[k] < 0.05),
      sig_analytic_same_sign = sum(same[k] & pa[k] < 0.05),
      in_top_elbow_ft = sum(rk[k] <= ELBOW & same[k]),
      median_rank_ft = median(rk[k]),
      median_absNES_main = median(abs(tab$NES_main[k])),
      median_absNES_ft = median(abs(nes[k])),
      median_NES_ratio = median(nes[k] / tab$NES_main[k]))
  }
  a <- rows[[length(rows) - 1]]; r <- rows[[length(rows)]]
  f1 <- fisher.test(matrix(c(a$same_sign, a$n - a$same_sign, r$same_sign, r$n - r$same_sign), 2))
  f2 <- fisher.test(matrix(c(a$sig_perm_same_sign, a$n - a$sig_perm_same_sign, r$sig_perm_same_sign, r$n - r$sig_perm_same_sign), 2))
  f3 <- fisher.test(matrix(c(a$in_top_elbow_ft, a$n - a$in_top_elbow_ft, r$in_top_elbow_ft, r$n - r$in_top_elbow_ft), 2))
  cat(sprintf("   Fisher activated vs repressed: same sign p = %.2g; significant and same sign p = %.2g; stays in the top %d p = %.2g\n",
              f1$p.value, f2$p.value, ELBOW, f3$p.value))
}
ar <- do.call(rbind, rows)
write.table(ar, file.path(OUT_DIR, "elbow_activated_vs_repressed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n== Elbow (%d): activated vs repressed ==\n", ELBOW))
print(ar %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)

cat("\n== Repressed regulators above the elbow that change sign or lose significance ==\n")
for (nm in main) {
  nes <- tab[[paste0("NES_", nm)]]; pe <- tab[[paste0("padj_emp_", nm)]]
  k <- tab$above_elbow & tab$direction_main == "repressed"
  cat(nm, "- change sign:", paste(sprintf("%s (%.1f -> %.1f)", tab$TF[k & nes > 0], tab$NES_main[k & nes > 0], nes[k & nes > 0]), collapse = ", "), "\n")
  cat(nm, "- same sign, not significant:", paste(tab$TF[k & nes < 0 & pe >= 0.05], collapse = ", "), "\n")
}
k <- tab$above_elbow & tab$direction_main == "repressed"
both <- k & tab$NES_GSE10971 < 0 & tab$NES_GSE69428 < 0 & tab$padj_emp_GSE10971 < 0.05 & tab$padj_emp_GSE69428 < 0.05
cat("Repressed, significant and with the same sign in BOTH series:", sum(both), "of", sum(k), "\n")
k2 <- tab$above_elbow & tab$direction_main == "activated"
both2 <- k2 & tab$NES_GSE10971 > 0 & tab$NES_GSE69428 > 0 & tab$padj_emp_GSE10971 < 0.05 & tab$padj_emp_GSE69428 < 0.05
cat("Activated, significant and with the same sign in BOTH series:", sum(both2), "of", sum(k2), "\n")

cat(sprintf("\n== Overlap of the top %d by |NES| (main vs fallopian tube, shared regulators) ==\n", ELBOW))
for (nm in main) {
  rk <- tab[[paste0("rank_abs_NES_", nm)]]
  cat(sprintf("%s: %d of the %d above the elbow are in the top %d of the fallopian tube analysis\n", nm,
              sum(tab$above_elbow & rk <= ELBOW, na.rm = TRUE), ELBOW, ELBOW))
}

cat("\n== Most repressed regulators with the fallopian tube signature (Stouffer) and their main NES ==\n")
print(tab %>% arrange(NES_Stouffer) %>%
        select(TF, NES_main, rank_abs_NES, above_elbow, NES_GSE10971, NES_GSE69428, NES_Stouffer) %>%
        head(20) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)
cat("\n== Most activated regulators with the fallopian tube signature (Stouffer) ==\n")
print(tab %>% arrange(desc(NES_Stouffer)) %>%
        select(TF, NES_main, rank_abs_NES, above_elbow, NES_GSE10971, NES_GSE69428, NES_Stouffer) %>%
        head(10) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)

core <- read.delim(file.path(BASE, "5_MRA/signature_core.tsv"))
cat(sprintf("\n== Signature of the main MRA (%s genes): sign against fallopian tube ==\n", format(nrow(core), big.mark = ",")))
de <- readRDS(file.path(OUT_DIR, "de_all.rds")); st <- read.delim(file.path(OUT_DIR, "de_Stouffer.tsv"))
ef <- list(GSE10971 = setNames(de$GSE10971$t, de$GSE10971$gene), GSE69428 = setNames(de$GSE69428$t, de$GSE69428$gene),
           Stouffer = setNames(st$z, st$gene))
for (nm in main) {
  e <- ef[[nm]][core$gene]; ok <- !is.na(e)
  for (s in c(1, -1)) {
    k <- ok & sign(core$z_mean) == s
    cat(sprintf("%s, %s genes: %d measured, %d (%.1f%%) same sign; Spearman rho with z_mean (all) = %.3f\n", nm,
                ifelse(s > 0, "up", "down"), sum(k), sum(sign(e[k]) == s), 100 * mean(sign(e[k]) == s),
                cor(e[ok], core$z_mean[ok], method = "spearman")))
  }
}
write.table(tab, file.path(OUT_DIR, "ft_mra_vs_main_ranked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
