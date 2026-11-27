#!/usr/bin/env Rscript

# ============================================================
# Fallopian tube contrast, step 3: concordance of the consensus signature
# with the tumor vs fallopian tube contrasts.
#
# Consensus = rows of 3_0_2_allgenes.tsv with n_sources == 3 (unique
# genes), direction = consensus_direction. For each contrast: fraction of
# consensus genes with the same sign, binomial test, and an empirical
# null of random genes among the measured ones.
#
# Input:  3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv,
#         9_Revision_analyses/fallopian_tube/de_all.rds, de_Stouffer.tsv (9_4_2)
# Output (9_Revision_analyses/fallopian_tube/):
#         consensus_concordance.tsv, consensus_concordance_by_gene.tsv
# ============================================================

suppressPackageStartupMessages({library(dplyr)})
set.seed(20261126)

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/fallopian_tube")
N_NULL  <- 10000

allg <- read.delim(file.path(BASE, "3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv"))
cons <- allg %>% filter(n_sources == 3) %>% distinct(gene, consensus_direction)
stopifnot(!anyDuplicated(cons$gene))
cat("Consensus:", nrow(cons), "genes;", sum(cons$consensus_direction == "up"), "up,",
    sum(cons$consensus_direction == "down"), "down\n")

de <- readRDS(file.path(OUT_DIR, "de_all.rds"))
st <- read.delim(file.path(OUT_DIR, "de_Stouffer.tsv"))
# uniform table per contrast: gene, effect (sign), padj
contrasts <- lapply(de, function(d) data.frame(gene = d$gene, effect = d$logFC, padj = d$adj.P.Val))
contrasts$Stouffer <- data.frame(gene = st$gene, effect = st$z, padj = st$adj.P.gene)
main <- c("GSE10971", "GSE69428", "Stouffer")
contrasts <- contrasts[c(main, setdiff(names(contrasts), main))]

# concordance of a gene set with its expected direction (+1/-1)
conc <- function(effect, padj, expected) {
  ok <- sign(effect) == expected
  c(n = length(ok), same = sum(ok), sig_conc = sum(ok & padj < 0.05), sig_disc = sum(!ok & padj < 0.05))
}

rows <- list(); by_gene <- cons
for (nm in names(contrasts)) {
  d <- contrasts[[nm]]
  m <- inner_join(cons, d, by = "gene") %>% mutate(expected = ifelse(consensus_direction == "up", 1, -1))
  by_gene[[paste0("effect_", nm)]] <- d$effect[match(by_gene$gene, d$gene)]
  by_gene[[paste0("padj_", nm)]]   <- d$padj[match(by_gene$gene, d$gene)]
  for (dir in c("all", "up", "down")) {
    mm <- if (dir == "all") m else m[m$consensus_direction == dir, ]
    o  <- conc(mm$effect, mm$padj, mm$expected)
    bt <- binom.test(o["same"], o["n"], 0.5)
    # empirical null: the same number of random genes among the measured ones, with the same expected directions
    null <- replicate(N_NULL, {
      i  <- sample(nrow(d), nrow(mm))
      ok <- sign(d$effect[i]) == mm$expected
      c(mean(ok), mean(ok & d$padj[i] < 0.05))
    })
    rows[[length(rows) + 1]] <- data.frame(
      contrast = nm, consensus_direction = dir,
      n_consensus = if (dir == "all") nrow(cons) else sum(cons$consensus_direction == dir),
      n_measured = o[["n"]], n_same_sign = o[["same"]], frac_same_sign = o[["same"]] / o[["n"]],
      ci95_low = bt$conf.int[1], ci95_high = bt$conf.int[2], p_binomial = bt$p.value,
      null_frac_mean = mean(null[1, ]), null_frac_p2.5 = quantile(null[1, ], 0.025),
      null_frac_p97.5 = quantile(null[1, ], 0.975),
      p_empirical = (1 + sum(null[1, ] >= o[["same"]] / o[["n"]])) / (N_NULL + 1),
      n_sig_concordant = o[["sig_conc"]], frac_sig_concordant = o[["sig_conc"]] / o[["n"]],
      null_frac_sig_mean = mean(null[2, ]),
      p_empirical_sig = (1 + sum(null[2, ] >= o[["sig_conc"]] / o[["n"]])) / (N_NULL + 1),
      n_sig_discordant = o[["sig_disc"]], row.names = NULL)
  }
}
tab <- do.call(rbind, rows)
write.table(tab, file.path(OUT_DIR, "consensus_concordance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(by_gene, file.path(OUT_DIR, "consensus_concordance_by_gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250)
print(tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# concordance in the two main series at the same time
b <- by_gene %>% filter(!is.na(effect_GSE10971), !is.na(effect_GSE69428)) %>%
  mutate(expected = ifelse(consensus_direction == "up", 1, -1),
         both = sign(effect_GSE10971) == expected & sign(effect_GSE69428) == expected,
         neither = sign(effect_GSE10971) != expected & sign(effect_GSE69428) != expected)
cat("\nConsensus genes measured in both series:", nrow(b), "\n")
print(b %>% group_by(consensus_direction) %>%
        summarise(n = n(), concordant_in_both = sum(both), discordant_in_both = sum(neither)))

# gradient: genes significant in 1, 2 or 3 sources of the main analysis (consensus direction)
cat("\nConcordance by number of sources of the main analysis (Stouffer and each series):\n")
ag <- allg %>% distinct(gene, n_sources, consensus_direction) %>% filter(consensus_direction != "none")
for (nm in main) {
  d <- contrasts[[nm]]
  g <- inner_join(ag, d, by = "gene") %>%
    mutate(ok = sign(effect) == ifelse(consensus_direction == "up", 1, -1)) %>%
    group_by(n_sources, consensus_direction) %>%
    summarise(n = n(), frac_same_sign = round(mean(ok), 3), .groups = "drop")
  cat("--", nm, "\n"); print(as.data.frame(g))
}
# overall fraction of genes with a positive effect (reference for the null)
for (nm in main) cat(nm, ": fraction of measured genes with effect > 0 =", round(mean(contrasts[[nm]]$effect > 0), 3), "\n")
sessionInfo()
