#!/usr/bin/env Rscript

# ============================================================
# Overall survival associated with the candidate transcriptional drivers.
#
# The drivers are the master regulators that are enriched in the
# highly dependent fraction of the consensus signature and are
# themselves essential in ovarian cancer cell lines (median CRISPR
# dependency <= -0.5), i.e. the same set drawn in Figure 4.
#
# Every driver is tested on its own, never as a combined score:
#   - Cox proportional hazards on the expression as a continuous
#     variable (hazard ratio per standard deviation).
#   - Cox adjusted for age, clinical stage and residual disease.
#   - Cox adjusted for the stromal and immune content of the tumour,
#     the two microenvironment axes that dominate the prognosis of
#     high-grade serous carcinoma and are the obvious confounder of
#     any bulk expression measurement.
#   - Median split of the cohort with a log-rank test.
#
# A permutation test that keeps the correlation between the drivers
# intact reports whether the candidate behave as a group or whether the
# signal is gene specific.
#
# The forest plot reports all of them; the Kaplan-Meier curves are
# drawn only for the drivers that survive the multiple-testing
# correction.
#
# Outputs: 8_1_0_cox_drivers_survival.tsv
#          8_1_1_forest_cox_drivers.pdf / .png
#          8_1_2_KM_significant_drivers.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
out_dir  <- file.path(base_dir, "8_Survival_analysis")

paths <- list(
  clinical  = file.path(out_dir, "8_0_clinical_TCGA_OV.rds"),
  counts    = file.path(base_dir, "3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_2_norm_counts_tumor_ONLY.rds"),
  metadata  = file.path(base_dir, "3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_1_metadata_all.rds"),
  annot     = file.path(base_dir, "3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_2_gene_annotation_UNION_universe.rds"),
  dep_cache = file.path(base_dir, "7_Figures/7_2_depmap_ovary_24Q4_646_cache.rds"),
  mra_cache = file.path(base_dir, "7_Figures/7_2_mra_enrichment_24Q4_646_cache.rds")
)

DEP_MR  <- -0.5    # essentiality threshold that defines the drivers
MIN_TPM <- 1       # median expression required in the same cell lines
P_SIG   <- 0.05

BASE_SIZE <- 15
col_up   <- "#B2182B"   # higher expression, shorter survival
col_down <- "#2166AC"   # higher expression, longer survival
col_ns   <- "#9E9E9E"

read_rds_here <- function(p) readRDS(path.expand(p))

# ============================================================
# 1) The candidate transcriptional drivers
# ============================================================
dep <- read_rds_here(paths$dep_cache)
mra <- read_rds_here(paths$mra_cache)
# the enrichment cache of 7_2 must be newer than the regulator table it derives from
mra_tsv <- path.expand(file.path(base_dir, "5_MRA/core_mra.tsv"))
stopifnot("rerun 7_2: its enrichment cache is older than 5_MRA/core_mra.tsv" =
            file.mtime(path.expand(paths$mra_cache)) > file.mtime(mra_tsv))
stopifnot("rerun 7_2: the dependency cache has no expression data" = !is.null(dep$expr_tpm))

drivers <- dep$crispr %>%
  filter(gene %in% mra$results$regulon) %>%
  group_by(gene) %>%
  summarise(median_dep = median(dep, na.rm = TRUE), .groups = "drop") %>%
  filter(median_dep <= DEP_MR) %>%
  filter(gene %in% names(dep$expr_tpm)[dep$expr_tpm >= MIN_TPM]) %>%
  arrange(median_dep)

cat("[INFO] transcriptional drivers tested:", nrow(drivers), "\n")

# ============================================================
# 2) Expression at the patient level
# ============================================================
counts   <- read_rds_here(paths$counts)
metadata <- read_rds_here(paths$metadata)
annot    <- read_rds_here(paths$annot)

# Primary tumours only: recurrences would mix different points of
# the disease course into a single overall-survival model
primary_ids <- metadata %>%
  filter(sample.type == "primary", sample.id %in% colnames(counts)) %>%
  pull(sample.id)

# The annotation of the count matrix predates a few HGNC renames
gene_aliases <- c(H2AZ1 = "H2AFZ", PHB1 = "PHB")

driver_lookup <- drivers %>%
  transmute(
    gene,
    annot_symbol = ifelse(gene %in% names(gene_aliases), gene_aliases[gene], gene)
  )

sym2ens <- annot %>%
  dplyr::select(ensembl, annot_symbol = gene) %>%
  inner_join(driver_lookup, by = "annot_symbol") %>%
  distinct(gene, .keep_all = TRUE)

expr <- counts[rownames(counts) %in% sym2ens$ensembl, primary_ids, drop = FALSE]
rownames(expr) <- sym2ens$gene[match(rownames(expr), sym2ens$ensembl)]

# Microenvironment scores: mean of the per-gene z-scores of a stromal
# and an immune marker set, summarised per patient the same way as the
# drivers. They are the covariates of the sensitivity analysis.
signature_score <- function(symbols) {
  ens <- intersect(annot$ensembl[annot$gene %in% symbols], rownames(counts))
  m   <- log2(counts[ens, primary_ids, drop = FALSE] + 1)
  v   <- colMeans(t(scale(t(m))), na.rm = TRUE)
  tapply(v, substr(primary_ids, 1, 12), median)
}

stromal_raw <- signature_score(c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "FBN1",
                                 "DCN", "LUM", "THBS2", "POSTN", "FAP", "ZFHX4",
                                 "SPARC", "VCAN", "MMP2"))
immune_raw  <- signature_score(c("PTPRC", "CD3D", "CD3E", "CD2", "CXCL9",
                                 "CXCL10", "TAP1", "GZMB", "IL2RG", "LCK"))

genes_ok <- intersect(drivers$gene, rownames(expr))
if (length(setdiff(drivers$gene, genes_ok)) > 0)
  stop("Drivers missing from the expression matrix: ",
       paste(setdiff(drivers$gene, genes_ok), collapse = ", "))
cat("[INFO] drivers found in the expression matrix:", length(genes_ok),
    "| primary tumours:", ncol(expr), "\n")

# Counts are library-size normalised; the log keeps the hazard ratio
# on an interpretable scale and limits the leverage of extreme values
expr_patient <- log2(expr[genes_ok, , drop = FALSE] + 1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample.id") %>%
  mutate(patient = substr(sample.id, 1, 12)) %>%
  group_by(patient) %>%
  summarise(across(all_of(genes_ok), ~ median(.x, na.rm = TRUE)), .groups = "drop")

# ============================================================
# 3) Analysis set
# ============================================================
clinical <- read_rds_here(paths$clinical)

dat <- expr_patient %>%
  inner_join(clinical, by = c("patient" = "bcr_patient_barcode")) %>%
  filter(!is.na(overall_survival), !is.na(deceased)) %>%
  mutate(
    os_months      = overall_survival / 30.44,
    stromal_score  = as.numeric(scale(stromal_raw[patient])),
    immune_score   = as.numeric(scale(immune_raw[patient])),
    stage_simple   = droplevels(factor(as.character(stage_simple))),
    residual_group = droplevels(factor(as.character(residual_group)))
  )

# Stage I and II are too small to be their own strata in this cohort
dat <- dat %>%
  mutate(stage_grouped = factor(
    ifelse(stage_simple %in% c("Stage I", "Stage II"), "Stage I/II",
           as.character(stage_simple)),
    levels = c("Stage I/II", "Stage III", "Stage IV")
  ))

cat("[INFO] patients analysed:", nrow(dat),
    "| deaths:", sum(dat$deceased), "\n")

surv_obj <- Surv(dat$os_months, dat$deceased)

# Median follow-up needs the reverse Kaplan-Meier estimator: the plain
# median of the observed times mixes deaths with censoring and badly
# underestimates how long the cohort was actually followed
rev_km <- survfit(Surv(dat$os_months, !dat$deceased) ~ 1)
cat("[INFO] median overall survival (months):",
    round(unname(summary(survfit(surv_obj ~ 1))$table["median"]), 1), "\n")
cat("[INFO] median follow-up, reverse Kaplan-Meier (months):",
    round(unname(summary(rev_km)$table["median"]), 1), "\n")

# ============================================================
# 4) One model per driver
# ============================================================
test_driver <- function(gene) {

  z    <- as.numeric(scale(dat[[gene]]))
  cut  <- median(dat[[gene]], na.rm = TRUE)
  grp  <- factor(ifelse(dat[[gene]] > cut, "High", "Low"), levels = c("Low", "High"))

  # --- continuous, unadjusted -------------------------------------
  fit  <- coxph(surv_obj ~ z)
  s    <- summary(fit)
  zph  <- tryCatch(cox.zph(fit)$table["z", "p"], error = function(e) NA_real_)

  # --- continuous, adjusted for the standard prognostic factors ---
  adj_df <- dat %>%
    mutate(z = z) %>%
    dplyr::select(os_months, deceased, z, age_dx, stage_grouped, residual_group)

  fit_adj <- tryCatch(
    coxph(Surv(os_months, deceased) ~ z + age_dx + stage_grouped + residual_group,
          data = adj_df),
    error = function(e) NULL
  )
  s_adj <- if (is.null(fit_adj)) NULL else summary(fit_adj)

  # --- continuous, adjusted for the tumour microenvironment ---------
  micro_df <- dat %>%
    mutate(z = z) %>%
    dplyr::select(os_months, deceased, z, stromal_score, immune_score)

  s_micro <- summary(coxph(Surv(os_months, deceased) ~ z + stromal_score + immune_score,
                           data = micro_df))

  # --- median split -----------------------------------------------
  fit_grp <- coxph(surv_obj ~ grp)
  s_grp   <- summary(fit_grp)
  lr      <- survdiff(surv_obj ~ grp)
  lr_p    <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

  med <- summary(survfit(surv_obj ~ grp))$table[, "median"]

  tibble(
    gene          = gene,
    n             = s$n,
    n_events      = s$nevent,
    HR_per_SD     = unname(s$conf.int[1, "exp(coef)"]),
    HR_lo95       = unname(s$conf.int[1, "lower .95"]),
    HR_hi95       = unname(s$conf.int[1, "upper .95"]),
    z_stat        = unname(s$coefficients[1, "z"]),
    p_value       = unname(s$coefficients[1, "Pr(>|z|)"]),
    ph_p          = zph,
    HR_adj        = if (is.null(s_adj)) NA_real_ else unname(s_adj$conf.int["z", "exp(coef)"]),
    HR_adj_lo95   = if (is.null(s_adj)) NA_real_ else unname(s_adj$conf.int["z", "lower .95"]),
    HR_adj_hi95   = if (is.null(s_adj)) NA_real_ else unname(s_adj$conf.int["z", "upper .95"]),
    p_adjusted_model = if (is.null(s_adj)) NA_real_ else unname(s_adj$coefficients["z", "Pr(>|z|)"]),
    HR_micro         = unname(s_micro$conf.int["z", "exp(coef)"]),
    HR_micro_lo95    = unname(s_micro$conf.int["z", "lower .95"]),
    HR_micro_hi95    = unname(s_micro$conf.int["z", "upper .95"]),
    p_micro          = unname(s_micro$coefficients["z", "Pr(>|z|)"]),
    HR_high_vs_low   = unname(s_grp$conf.int[1, "exp(coef)"]),
    HR_high_lo95     = unname(s_grp$conf.int[1, "lower .95"]),
    HR_high_hi95     = unname(s_grp$conf.int[1, "upper .95"]),
    logrank_p        = lr_p,
    median_OS_low    = unname(med[1]),
    median_OS_high   = unname(med[2])
  )
}

res <- map_dfr(genes_ok, test_driver) %>%
  mutate(
    p_bh          = p.adjust(p_value, method = "BH"),
    p_micro_bh    = p.adjust(p_micro, method = "BH"),
    logrank_p_bh  = p.adjust(logrank_p, method = "BH"),
    direction     = case_when(
      p_value >= P_SIG ~ "Not significant",
      HR_per_SD > 1    ~ "Higher expression, shorter survival",
      TRUE             ~ "Higher expression, longer survival"
    )
  ) %>%
  left_join(drivers, by = "gene") %>%
  arrange(p_value)

write_tsv(res, file.path(path.expand(out_dir), "8_1_0_cox_drivers_survival.tsv"))

cat("\n[RESULT] nominally significant (p < 0.05):",
    sum(res$p_value < P_SIG), "of", nrow(res), "\n")
cat("[RESULT] significant after BH correction:", sum(res$p_bh < P_SIG), "\n\n")
print(as.data.frame(res %>%
  dplyr::select(gene, HR_per_SD, HR_lo95, HR_hi95, p_value, p_bh,
                logrank_p, HR_adj, p_adjusted_model, HR_micro, p_micro_bh) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# ============================================================
# 4b) Do the candidate behave as a group?
# ============================================================
# The drivers are co-expressed, so their hazard ratios are not
# independent observations. Permuting the survival labels keeps that
# correlation intact and gives the honest null for the question of
# whether the set as a whole leans towards protection.
N_PERM <- 1000
set.seed(11)

z_mat <- vapply(genes_ok, function(g) as.numeric(scale(dat[[g]])),
                numeric(nrow(dat)))

frac_protective <- function(time, event) {
  suppressWarnings(
    mean(vapply(seq_len(ncol(z_mat)),
                function(j) coef(coxph(Surv(time, event) ~ z_mat[, j]))[1],
                numeric(1)) < 0)
  )
}

obs_frac  <- frac_protective(dat$os_months, dat$deceased)
null_frac <- replicate(N_PERM, {
  i <- sample(nrow(dat))
  frac_protective(dat$os_months[i], dat$deceased[i])
})

cat("\n[GROUP] drivers with HR < 1:", round(obs_frac, 3),
    "| permutation null:", round(mean(null_frac), 3),
    "| p =", round(mean(null_frac >= obs_frac), 4), "\n\n")

# ============================================================
# 5) Forest plot
# ============================================================
pal <- c(
  "Higher expression, shorter survival" = col_up,
  "Higher expression, longer survival"  = col_down,
  "Not significant"                     = col_ns
)

# Only the statistic the conclusion rests on is printed: the hazard
# ratio and the p value after the correction across the candidate tests
forest_df <- res %>%
  mutate(
    gene  = factor(gene, levels = rev(res$gene[order(res$HR_per_SD)])),
    label = sprintf("%.2f (%.2f-%.2f)", HR_per_SD, HR_lo95, HR_hi95),
    p_lab = ifelse(p_bh < 0.001, "< 0.001", sprintf("%.3f", p_bh))
  )

# Both halves of the forest share this expansion so that their rows
# stay aligned once patchwork puts them side by side
y_expand <- expansion(add = c(0.6, 1.5))

p_forest <- ggplot(forest_df, aes(x = HR_per_SD, y = gene, colour = direction)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(xmin = HR_lo95, xmax = HR_hi95), orientation = "y",
                width = 0.28, linewidth = 0.7) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_y_discrete(limits = levels(forest_df$gene), expand = y_expand) +
  scale_x_continuous(trans = "log2", breaks = c(0.7, 0.85, 1, 1.2, 1.4)) +
  labs(x = "Hazard ratio per standard deviation of expression", y = NULL) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    axis.text.y     = element_text(face = "bold", colour = "grey10"),
    axis.title.x    = element_text(face = "bold"),
    legend.position = "top",
    legend.text     = element_text(size = BASE_SIZE - 3),
    plot.margin     = margin(12, 6, 12, 12)
  )

# Numbers printed next to the forest, on their own strip
p_table <- ggplot(forest_df, aes(y = gene)) +
  geom_text(aes(x = 0, label = label), hjust = 0, size = 4.1, colour = "grey15") +
  geom_text(aes(x = 1, label = p_lab), hjust = 0, size = 4.1, colour = "grey15") +
  annotate("text", x = 0, y = nrow(forest_df) + 0.9, label = "HR (95% CI)",
           hjust = 0, size = 4.3, fontface = "bold") +
  annotate("text", x = 1, y = nrow(forest_df) + 0.9, label = "BH-adjusted p",
           hjust = 0, size = 4.3, fontface = "bold") +
  scale_x_continuous(limits = c(0, 1.75)) +
  scale_y_discrete(limits = levels(forest_df$gene), expand = y_expand) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = BASE_SIZE) +
  theme(plot.margin = margin(12, 12, 12, 0))

figure_forest <- p_forest + p_table + plot_layout(widths = c(2.6, 1))

ggsave(file.path(path.expand(out_dir), "8_1_1_forest_cox_drivers.pdf"),
       figure_forest, width = 11.5, height = 7.2, device = cairo_pdf)
ggsave(file.path(path.expand(out_dir), "8_1_1_forest_cox_drivers.png"),
       figure_forest, width = 11.5, height = 7.2, dpi = 400, bg = "white")

# ============================================================
# 6) Kaplan-Meier curves of the drivers that survive the correction
# ============================================================
km_genes <- res %>% filter(p_bh < P_SIG) %>% arrange(p_value) %>% pull(gene)

cat("[INFO] Kaplan-Meier curves drawn for:", paste(km_genes, collapse = ", "), "\n")

# The axis has to cover the whole follow-up: truncating it earlier
# hides the flat tail of a censored patient while the numbers at risk,
# which come from the fit and not from the axis, still count them
MAX_TIME   <- max(dat$os_months, na.rm = TRUE)
RISK_TIMES <- seq(0, floor(MAX_TIME / 30) * 30, by = 30)
X_LIM      <- c(0, MAX_TIME * 1.02)

km_panel <- function(g) {

  cut <- median(dat[[g]], na.rm = TRUE)
  grp <- factor(ifelse(dat[[g]] > cut, "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(dat$os_months, dat$deceased) ~ grp)

  strata <- rep(names(fit$strata), fit$strata)
  curve <- tibble(
    group    = factor(sub("grp=", "", strata), levels = c("Low", "High")),
    time     = fit$time,
    surv     = fit$surv,
    n_censor = fit$n.censor
  ) %>%
    group_by(group) %>%
    # the curve has to start at 100 % survival on day zero
    group_modify(~ bind_rows(tibble(time = 0, surv = 1, n_censor = 0), .x)) %>%
    ungroup()

  st <- res %>% filter(gene == g)

  # Each raw p is shown next to the correction of its own family: the
  # log-rank belongs to the median split drawn here, whereas the BH
  # value quoted in the text comes from the continuous Cox models,
  # which is the test the candidate drivers were screened with
  fmt_p <- function(p, prefix = "p") {
    ifelse(p < 0.001, paste(prefix, "< 0.001"), sprintf("%s = %.3f", prefix, p))
  }

  lab <- sprintf(
    paste0("HR = %.2f per SD, 95%% CI %.2f-%.2f, BH-adjusted %s\n",
           "Median split shown for display, log-rank %s"),
    st$HR_per_SD, st$HR_lo95, st$HR_hi95, fmt_p(st$p_bh),
    fmt_p(st$logrank_p)
  )

  p_curve <- ggplot(curve, aes(time, surv, colour = group)) +
    geom_step(linewidth = 0.9) +
    geom_point(data = filter(curve, n_censor > 0), shape = 3, size = 1.6,
               show.legend = FALSE) +
    scale_colour_manual(values = c(Low = col_down, High = col_up),
                        name = "Expression",
                        labels = c("Low (below median)", "High (above median)")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_x_continuous(limits = X_LIM, breaks = RISK_TIMES,
                       expand = expansion(mult = 0.07)) +
    # The statistics go above the panel: three lines of text do not fit
    # inside it without running into one of the two curves
    labs(title = g, subtitle = lab, x = NULL, y = "Overall survival probability") +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title       = element_text(face = "bold", size = BASE_SIZE + 3, hjust = 0),
      plot.subtitle    = element_text(size = BASE_SIZE - 3, colour = "grey25",
                                      lineheight = 1.2, margin = margin(4, 0, 8, 0)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
      axis.title       = element_text(face = "bold"),
      legend.title     = element_text(face = "bold"),
      plot.margin      = margin(10, 14, 4, 10)
    )

  # Numbers at risk with the cumulative deaths accumulated up to each
  # time point, printed under the corresponding curve
  s <- summary(fit, times = RISK_TIMES, extend = TRUE)

  risk <- tibble(
    group   = factor(sub("grp=", "", s$strata), levels = c("Low", "High")),
    time    = s$time,
    n_risk  = s$n.risk,
    n_event = s$n.event
  ) %>%
    group_by(group) %>%
    mutate(cum_event = cumsum(n_event)) %>%
    ungroup()

  number_track <- function(df, value, title, axis_title = NULL) {
    ggplot(df, aes(time, group, colour = group)) +
      geom_text(aes(label = {{ value }}), size = 3.9, show.legend = FALSE) +
      scale_colour_manual(values = c(Low = col_down, High = col_up)) +
      scale_x_continuous(limits = X_LIM, breaks = RISK_TIMES,
                         expand = expansion(mult = 0.07)) +
      scale_y_discrete(limits = c("High", "Low")) +
      labs(title = title, x = axis_title, y = NULL) +
      theme_bw(base_size = BASE_SIZE) +
      theme(
        plot.title   = element_text(size = BASE_SIZE - 2, colour = "grey25"),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        axis.ticks   = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_text(face = "bold", colour = "grey15"),
        axis.title.x = element_text(face = "bold"),
        plot.margin  = margin(2, 14, 2, 10)
      )
  }

  p_risk <- number_track(risk, n_risk, "Number at risk",
                         axis_title = "Overall survival (months)") +
    theme(plot.margin = margin(2, 14, 10, 10))

  p_curve / p_risk + plot_layout(heights = c(4.2, 1))
}

# Kaplan-Meier figure only when at least one driver passes the BH threshold
if (length(km_genes) == 0) {
  cat("[INFO] no driver is significant after BH correction; Kaplan-Meier figure not drawn\n")
  unlink(file.path(path.expand(out_dir), paste0("8_1_2_KM_significant_drivers.", c("pdf", "png"))))
} else {
km_panels <- map(km_genes, km_panel)

figure_km <- wrap_plots(km_panels, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

km_w <- 6.6 * length(km_genes)

ggsave(file.path(path.expand(out_dir), "8_1_2_KM_significant_drivers.pdf"),
       figure_km, width = km_w, height = 7.2, device = cairo_pdf)
ggsave(file.path(path.expand(out_dir), "8_1_2_KM_significant_drivers.png"),
       figure_km, width = km_w, height = 7.2, dpi = 400, bg = "white")

message("Survival analysis saved to: ", out_dir)
}
