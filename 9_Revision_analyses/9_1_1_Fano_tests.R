#!/usr/bin/env Rscript

# ============================================================
# Statistical tests of the Fano factor by gene subset.
#
# Compares the consensus genes with the other subsets of differentially
# expressed genes, on the Fano factor and on its residual after removing
# the mean-variance trend, in each control cohort (GTEx and AE). Also
# tests the difference between cohorts, paired by gene, within each subset.
#
# Input:  3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_3_fano_{GTEX,AE}_raw.rds
#         3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv
# Output: 9_Revision_analyses/consensus/fano_tests.tsv
#         9_Revision_analyses/consensus/fano_paired_cohorts.tsv
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(splines)
})

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
IN_DIR  <- file.path(BASE, "3_Consensus_DGE_analysis")
OUT_DIR <- file.path(BASE, "9_Revision_analyses/consensus")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Cliff's delta from the Wilcoxon statistic
cliff <- function(x, y) 2 * unname(wilcox.test(x, y)$statistic) / (length(x) * length(y)) - 1

fano <- list(
  GTEx = readRDS(file.path(IN_DIR, "3_3_Fano_Factor/3_3_3_fano_GTEX_raw.rds")),
  AE   = readRDS(file.path(IN_DIR, "3_3_Fano_Factor/3_3_3_fano_AE_raw.rds"))
)

all_genes <- read_tsv(file.path(IN_DIR, "3_1_Core_signature/3_0_2_allgenes.tsv"), show_col_types = FALSE) %>%
  filter(!is.na(ensembl)) %>%
  distinct(ensembl, .keep_all = TRUE) %>%
  mutate(subset = case_when(
    n_sources == 3        ~ "consensus",
    gtex == 1 & ae == 1   ~ "GTEx_AE",
    ae == 1 & geo == 1    ~ "AE_GEO",
    gtex == 1 & geo == 1  ~ "GTEx_GEO",
    TRUE                  ~ "one_source"
  ))

res <- list()
out <- list()
for (co in names(fano)) {
  d <- fano[[co]] %>%
    left_join(all_genes %>% select(gene = ensembl, subset), by = "gene") %>%
    mutate(subset = coalesce(subset, "not_DE"),
           lf  = log10(fano_factor + 1),
           lm_ = log10(mean_expression + 1))
  # Residual of the Fano factor on the mean (mean-variance trend of counts)
  d$res <- residuals(lm(lf ~ ns(lm_, 5), data = d))
  out[[co]] <- d %>% mutate(cohort = co)

  cat("\n==", co, "genes:", nrow(d), "\n")
  print(d %>% group_by(subset) %>%
          summarise(n = n(), med_fano = median(fano_factor), med_mean = median(mean_expression),
                    med_res = median(res)) %>%
          arrange(desc(n)))

  x <- d[d$subset == "consensus", ]
  for (ref in c("GTEx_AE", "one_source", "not_DE", "all_other_DE")) {
    y <- if (ref == "all_other_DE") d[!d$subset %in% c("consensus", "not_DE"), ] else d[d$subset == ref, ]
    res[[length(res) + 1]] <- tibble(
      cohort = co, comparison = paste("consensus vs", ref), n_cons = nrow(x), n_ref = nrow(y),
      med_fano_cons = median(x$fano_factor), med_fano_ref = median(y$fano_factor),
      p_fano = wilcox.test(x$fano_factor, y$fano_factor)$p.value,
      cliff_fano = cliff(x$fano_factor, y$fano_factor),
      p_residual = wilcox.test(x$res, y$res)$p.value,
      cliff_residual = cliff(x$res, y$res))
  }
}
res <- bind_rows(res)
print(as.data.frame(res), digits = 3)
write_tsv(res, file.path(OUT_DIR, "fano_tests.tsv"))

# Difference between cohorts, paired by gene, within each subset
w <- inner_join(out$GTEx %>% select(gene, subset, f_gtex = fano_factor),
                out$AE %>% select(gene, f_ae = fano_factor), by = "gene")
paired <- w %>%
  group_by(subset) %>%
  summarise(n = n(), med_gtex = median(f_gtex), med_ae = median(f_ae),
            delta_med = median(f_gtex) - median(f_ae),
            med_paired_diff = median(f_gtex - f_ae),
            p_signed = wilcox.test(f_gtex, f_ae, paired = TRUE)$p.value)
print(as.data.frame(paired), digits = 4)
write_tsv(paired, file.path(OUT_DIR, "fano_paired_cohorts.tsv"))
