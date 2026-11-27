#!/usr/bin/env Rscript

# ============================================================
# Overall survival of the candidate regulators in cohorts independent of
# TCGA (curatedOvarianData): Cox model per gene and dataset, meta-analysis
# per gene and group tests. No signatures, no cut-points.
#
# Input:  6_Depmap_ovary/6_3_2_candidate_regulators.tsv,
#         8_Survival_analysis/8_1_0_cox_drivers_survival.tsv,
#         curatedOvarianData and msigdbr (packages)
# Output (9_Revision_analyses/survival_independent/):
#         datasets_used.tsv, cox_by_dataset.tsv, meta_by_gene.tsv,
#         meta_sensitivity.tsv, leave_one_out.tsv, group_tests.tsv,
#         comparison_TCGA.tsv, aliases_used.tsv, results.rds
# ============================================================

suppressPackageStartupMessages({
  library(curatedOvarianData)
  library(Biobase)
  library(survival)
  library(data.table)
  library(parallel)
})

options(stringsAsFactors = FALSE, width = 200)
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

# ===================== CONFIG ===================== #
BASE     <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
out_dir  <- file.path(BASE, "9_Revision_analyses/survival_independent")
tcga_tsv <- file.path(BASE, "8_Survival_analysis/8_1_0_cox_drivers_survival.tsv")
cand_tsv <- file.path(BASE, "6_Depmap_ovary/6_3_2_candidate_regulators.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_CORES         <- 10     # the random draws depend on the number of workers: keep it to reproduce the nulls
MIN_DEATHS      <- 40     # main inclusion criterion
MIN_DEATHS_LAX  <- 15     # lax threshold (the default of the package), sensitivity only
MIN_GENES       <- 1000   # drops targeted panels (not genome-wide)
MIN_NON_NA      <- 0.8    # minimum fraction of patients with an expression value
N_RAND          <- 10000  # random sets of genes
N_PERM          <- 10000  # permutations of the survival labels
SEED            <- 11

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

# Candidate regulators of the pipeline, in alphabetical order
genes <- sort(fread(cand_tsv)[candidate == TRUE]$gene, method = "radix")

# Aliases tried, in order, when the official symbol is missing.
# AIM1 is not used for AURKB: on these arrays it is the symbol of another gene (CRYBG1).
alias <- list(
  H2AZ1    = c("H2AFZ", "H2AZ"),
  MIS18BP1 = c("C14orf106", "KNL2", "M18BP1"),
  CDK1     = c("CDC2"),
  SRSF2    = c("SFRS2"),
  AURKB    = c("STK12", "AIK2", "ARK2"),
  CHEK1    = c("CHK1"),
  ILF2     = c("NF45"),
  TAF5     = c("TAF2D"),
  WDHD1    = c("AND-1", "CTF4"),
  ACTL6A   = c("BAF53A"),
  PLK1     = c("PLK"),
  CCNA2    = c("CCNA"),
  TTF2     = c("HuF2")
)

tcga_sets <- c("TCGA_eset", "TCGA.RNASeqV2_eset", "TCGA.mirna.8x15kv2_eset")

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ============================================================
# 1) Load the datasets and remove duplicates
# ============================================================
# Blacklist of samples duplicated between studies (Waldron et al. 2014),
# as distributed with the package
cfg <- new.env()
sys.source(system.file("extdata", "patientselection.config", package = "curatedOvarianData"),
           envir = cfg)
blacklist <- cfg$duplicates

all_sets <- data(package = "curatedOvarianData")$results[, "Item"]

esets <- list()
info  <- list()
for (nm in all_sets) {
  e <- new.env()
  data(list = nm, package = "curatedOvarianData", envir = e)
  x <- e[[nm]]
  warns  <- paste(unlist(experimentData(x)@other$warnings), collapse = " ")
  drop   <- sub(".*:", "", blacklist[sub(":.*", "", blacklist) == nm])
  n_orig <- ncol(x)
  x <- x[, !sampleNames(x) %in% drop]
  # only features that map to a single gene
  x <- x[!grepl("///", featureNames(x), fixed = TRUE), ]
  esets[[nm]] <- x
  info[[nm]]  <- data.table(
    dataset = nm, platform = annotation(x)[1], pmid = pubMedIds(x)[1],
    n_samples = n_orig, n_blacklisted = length(drop), n_unique_genes = nrow(x),
    package_warning = substr(warns, 1, 60)
  )
}
info <- rbindlist(info)

# ============================================================
# 2) Subset analyzed in each dataset
# ============================================================
# Returns time, event, covariates and the expression matrix (genes x patients)
prepare <- function(nm, subset = c("serous", "high_grade_2_3", "high_grade_package")) {
  subset <- match.arg(subset)
  x <- esets[[nm]]
  p <- pData(x)
  ok <- p$sample_type %in% "tumor" & p$histological_type %in% "ser" &
    !is.na(p$days_to_death) & !is.na(p$vital_status) & p$days_to_death > 0
  if (subset == "high_grade_2_3")
    ok <- ok & (p$grade %in% 2:4 | (is.na(p$grade) & p$summarygrade %in% "high"))
  if (subset == "high_grade_package")
    ok <- ok & p$summarygrade %in% "high"
  # one tumor per patient
  id <- p$unique_patient_ID
  ok <- ok & !(duplicated(id) & !is.na(id))
  p  <- p[ok, , drop = FALSE]
  X  <- exprs(x)[, ok, drop = FALSE]

  # rename the aliases of the candidates to the official symbol
  used <- character(0)
  for (g in genes) {
    if (g %in% rownames(X)) { used[g] <- g; next }
    al <- intersect(alias[[g]], rownames(X))
    if (length(al)) {
      rownames(X)[rownames(X) == al[1]] <- g
      used[g] <- al[1]
    }
  }

  # stage grouped as in the TCGA analysis; a level with < 5 patients
  # is merged with stage III
  stg <- ifelse(p$tumorstage %in% 1:2, "I/II",
                ifelse(p$tumorstage %in% 3, "III", ifelse(p$tumorstage %in% 4, "IV", NA)))
  tb <- table(stg)
  stg[stg %in% names(tb)[tb < 5]] <- "III"
  stg <- factor(stg, levels = c("I/II", "III", "IV"))

  list(dataset = nm, subset = subset,
       time  = p$days_to_death / 30.44,
       event = as.integer(p$vital_status == "deceased"),
       age   = p$age_at_initial_pathologic_diagnosis,
       stage = droplevels(stg),
       X = X, symbol_used = used)
}

# ============================================================
# 3) Cox models
# ============================================================
# Fast univariate Cox (same estimate as coxph, without the wrapper)
cox_fast <- function(x, time, event) {
  ok <- !is.na(x)
  if (mean(ok) < MIN_NON_NA || sum(event[ok]) < 10) return(c(NA, NA, sum(ok), sum(event[ok])))
  s <- sd(x[ok])
  if (!is.finite(s) || s == 0) return(c(NA, NA, sum(ok), sum(event[ok])))
  z <- (x[ok] - mean(x[ok])) / s
  f <- tryCatch(
    coxph.fit(x = matrix(z, ncol = 1), y = Surv(time[ok], event[ok]), strata = NULL,
              offset = NULL, init = NULL, control = coxph.control(), weights = NULL,
              method = "efron", rownames = NULL, resid = FALSE),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(f)) return(c(NA, NA, sum(ok), sum(event[ok])))
  c(unname(f$coefficients), sqrt(f$var[1, 1]), sum(ok), sum(event[ok]))
}

# Cox of every gene of a dataset (universe for the random gene sets)
cox_genome <- function(d) {
  r <- mclapply(seq_len(nrow(d$X)), function(i) cox_fast(d$X[i, ], d$time, d$event),
                mc.cores = N_CORES)
  r <- do.call(rbind, r)
  data.table(dataset = d$dataset, gene = rownames(d$X),
             logHR = r[, 1], se = r[, 2], n = r[, 3], deaths = r[, 4])
}

# Detailed Cox of one candidate: unadjusted and adjusted for age and stage
cox_candidate <- function(d, g) {
  x  <- d$X[g, ]
  ok <- !is.na(x)
  df <- data.table(time = d$time, event = d$event, age = d$age, stage = d$stage,
                   z = NA_real_)
  df$z[ok] <- (x[ok] - mean(x[ok])) / sd(x[ok])
  df <- df[ok]
  fit <- coxph(Surv(time, event) ~ z, data = df)
  s   <- summary(fit)
  zph <- tryCatch(cox.zph(fit)$table["z", "p"], error = function(e) NA_real_)

  # covariates available in this dataset
  has_age   <- mean(!is.na(df$age)) >= 0.8
  has_stage <- mean(!is.na(df$stage)) >= 0.8 && nlevels(droplevels(df$stage)) > 1
  covs <- c(if (has_age) "age", if (has_stage) "stage")
  adj <- list(HR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, b = NA_real_,
              se = NA_real_, n = NA_integer_, ev = NA_integer_)
  if (length(covs)) {
    fa <- tryCatch(coxph(as.formula(paste("Surv(time, event) ~ z +", paste(covs, collapse = " + "))),
                         data = df), error = function(e) NULL)
    if (!is.null(fa)) {
      sa  <- summary(fa)
      adj <- list(HR = sa$conf.int["z", "exp(coef)"], lo = sa$conf.int["z", "lower .95"],
                  hi = sa$conf.int["z", "upper .95"], p = sa$coefficients["z", "Pr(>|z|)"],
                  b = sa$coefficients["z", "coef"], se = sa$coefficients["z", "se(coef)"],
                  n = sa$n, ev = sa$nevent)
    }
  }
  data.table(
    dataset = d$dataset, gene = g, feature_used = d$symbol_used[g],
    n = s$n, deaths = s$nevent,
    logHR = s$coefficients[1, "coef"], se = s$coefficients[1, "se(coef)"],
    HR = s$conf.int[1, "exp(coef)"], HR_lo95 = s$conf.int[1, "lower .95"],
    HR_hi95 = s$conf.int[1, "upper .95"], p = s$coefficients[1, "Pr(>|z|)"], p_PH = zph,
    adjustment_covariates = if (length(covs)) paste(covs, collapse = "+") else NA_character_,
    n_adj = adj$n, deaths_adj = adj$ev, logHR_adj = adj$b, se_adj = adj$se,
    HR_adj = adj$HR, HR_adj_lo95 = adj$lo, HR_adj_hi95 = adj$hi, p_adj = adj$p
  )
}

cox_candidates <- function(d) {
  g_ok <- intersect(genes, rownames(d$X))
  g_ok <- g_ok[vapply(g_ok, function(g) mean(!is.na(d$X[g, ])) >= MIN_NON_NA, logical(1))]
  rbindlist(lapply(g_ok, function(g) cox_candidate(d, g)))
}

# ============================================================
# 4) Meta-analysis: inverse variance (fixed) and DerSimonian-Laird (random)
# ============================================================
# Vectorized by rows: B and S are genes x datasets matrices, NA where there is no estimate
meta_mat <- function(B, S) {
  W  <- 1 / S^2
  W[is.na(B)] <- NA
  k  <- rowSums(!is.na(W))
  sw <- rowSums(W, na.rm = TRUE)
  bf <- rowSums(W * B, na.rm = TRUE) / sw
  Q  <- rowSums(W * (B - bf)^2, na.rm = TRUE)
  C  <- sw - rowSums(W^2, na.rm = TRUE) / sw
  tau2 <- ifelse(k > 1, pmax(0, (Q - (k - 1)) / C), 0)
  Wr <- 1 / (S^2 + tau2)
  Wr[is.na(B)] <- NA
  swr <- rowSums(Wr, na.rm = TRUE)
  br  <- rowSums(Wr * B, na.rm = TRUE) / swr
  data.table(
    k = k,
    logHR_fixed = bf, se_fixed = sqrt(1 / sw), p_fixed = 2 * pnorm(-abs(bf * sqrt(sw))),
    logHR_random = br, se_random = sqrt(1 / swr), p_random = 2 * pnorm(-abs(br * sqrt(swr))),
    Q = Q, p_Q = ifelse(k > 1, pchisq(Q, k - 1, lower.tail = FALSE), NA_real_),
    I2 = ifelse(k > 1 & Q > 0, pmax(0, 100 * (Q - (k - 1)) / Q), 0), tau2 = tau2
  )
}

# From a long table (dataset, gene, logHR, se) to a meta-analysis table per gene
meta_long <- function(dt, col_b = "logHR", col_se = "se", col_n = "n", col_ev = "deaths") {
  dt <- dt[!is.na(get(col_b))]
  B  <- dcast(dt, gene ~ dataset, value.var = col_b)
  S  <- dcast(dt, gene ~ dataset, value.var = col_se)
  m  <- meta_mat(as.matrix(B[, -1]), as.matrix(S[, -1]))
  tot <- dt[, .(n_patients = sum(get(col_n)), n_deaths = sum(get(col_ev)),
                n_HR_below_1 = sum(get(col_b) < 0),
                datasets = paste(sub("_eset$", "", dataset), collapse = ";")), by = gene]
  m <- cbind(gene = B$gene, m)
  m <- merge(m, tot, by = "gene")
  m[, `:=`(
    HR_fixed = exp(logHR_fixed), HR_fixed_lo95 = exp(logHR_fixed - 1.96 * se_fixed),
    HR_fixed_hi95 = exp(logHR_fixed + 1.96 * se_fixed),
    HR_random = exp(logHR_random), HR_random_lo95 = exp(logHR_random - 1.96 * se_random),
    HR_random_hi95 = exp(logHR_random + 1.96 * se_random),
    BH_fixed = p.adjust(p_fixed, "BH"), BH_random = p.adjust(p_random, "BH"))]
  m[order(p_random)]
}

# ============================================================
# 5) Which datasets are included
# ============================================================
counts <- rbindlist(lapply(setdiff(all_sets, tcga_sets), function(nm) {
  p  <- pData(esets[[nm]])
  d  <- prepare(nm, "serous")
  d2 <- prepare(nm, "high_grade_2_3")
  d3 <- prepare(nm, "high_grade_package")
  data.table(dataset = nm,
             n_tumor = sum(p$sample_type %in% "tumor"),
             n_serous = sum(p$sample_type %in% "tumor" & p$histological_type %in% "ser"),
             n_serous_OS = length(d$time), deaths_serous = sum(d$event),
             n_high_grade_2_3 = length(d2$time), deaths_high_grade_2_3 = sum(d2$event),
             n_high_grade_package = length(d3$time), deaths_high_grade_package = sum(d3$event),
             candidates_measured = sum(genes %in% rownames(d$X)),
             age_available = mean(!is.na(d$age)) >= 0.8 && length(d$age) > 0,
             stage_available = mean(!is.na(d$stage)) >= 0.8 && length(d$stage) > 0)
}))
dsets <- merge(info, counts, by = "dataset", all.x = TRUE)

dsets[, c("status", "reason") := {
  retract <- grepl("retracted", package_warning)
  subst   <- grepl("subset", package_warning)
  st <- fifelse(dataset %in% tcga_sets, "excluded",
        fifelse(subst, "excluded",
        fifelse(n_serous_OS == 0, "excluded",
        fifelse(n_unique_genes < MIN_GENES, "excluded",
        fifelse(retract, "sensitivity",
        fifelse(deaths_serous >= MIN_DEATHS, "main",
        fifelse(deaths_serous >= MIN_DEATHS_LAX, "sensitivity", "excluded")))))))
  rs <- fifelse(dataset %in% tcga_sets, "TCGA: discovery cohort",
        fifelse(subst, "subset of PMID17290060 according to the package",
        fifelse(n_serous_OS == 0, "no primary serous tumors with overall survival",
        fifelse(n_unique_genes < MIN_GENES, "targeted panel, not genome-wide",
        fifelse(retract, "retracted article: sensitivity only",
        fifelse(deaths_serous >= MIN_DEATHS, sprintf("at least %d deaths among serous tumors", MIN_DEATHS),
        fifelse(deaths_serous >= MIN_DEATHS_LAX,
                sprintf("fewer than %d deaths among serous tumors: sensitivity only (lax threshold %d)",
                        MIN_DEATHS, MIN_DEATHS_LAX),
                sprintf("fewer than %d deaths among serous tumors", MIN_DEATHS_LAX))))))))
  list(st, rs)
}]
dsets[, package_warning := NULL]
setorder(dsets, status, -deaths_serous)

sets_main    <- dsets[status == "main", dataset]
sets_retract <- dsets[status == "sensitivity" & grepl("retracted", reason), dataset]
sets_lax     <- dsets[status == "sensitivity" & grepl("lax threshold", reason), dataset]

msg("[INFO] datasets in the main analysis (%d): %s", length(sets_main),
    paste(sub("_eset", "", sets_main), collapse = ", "))
msg("[INFO] sensitivity only: %s", paste(sub("_eset", "", c(sets_retract, sets_lax)), collapse = ", "))

# ============================================================
# 6) Cox per dataset
# ============================================================
data_main <- lapply(setNames(sets_main, sets_main), prepare, subset = "serous")

# 6a) candidates, main analysis and sensitivity analyses
cox_main <- rbindlist(lapply(data_main, cox_candidates))
cox_main[, analysis := "main"]

cox_extra <- rbindlist(lapply(c(sets_retract, sets_lax), function(nm)
  cox_candidates(prepare(nm, "serous"))))
cox_extra[, analysis := fifelse(dataset %in% sets_retract, "extra_retracted", "extra_threshold_15")]

cox_sub <- rbindlist(lapply(c("high_grade_2_3", "high_grade_package"), function(sc) {
  ok <- dsets[status %in% c("main") & get(paste0("deaths_", sc)) >= MIN_DEATHS, dataset]
  r  <- rbindlist(lapply(ok, function(nm) cox_candidates(prepare(nm, sc))))
  r[, analysis := sc]
  r
}))

cox_all <- rbindlist(list(cox_main, cox_sub, cox_extra), use.names = TRUE)
setcolorder(cox_all, c("analysis", "dataset", "gene"))

# ============================================================
# 7) Meta-analysis per gene
# ============================================================
meta_main <- meta_long(cox_main)

sens <- list(
  main               = meta_main,
  high_grade_2_3     = meta_long(cox_sub[analysis == "high_grade_2_3"]),
  high_grade_package = meta_long(cox_sub[analysis == "high_grade_package"]),
  adjusted_age_and_stage = meta_long(cox_main[adjustment_covariates == "age+stage"],
                                     "logHR_adj", "se_adj", "n_adj", "deaths_adj"),
  adjusted_available_covariates = meta_long(cox_main[!is.na(adjustment_covariates)],
                                            "logHR_adj", "se_adj", "n_adj", "deaths_adj"),
  with_retracted      = meta_long(rbind(cox_main, cox_extra[analysis == "extra_retracted"])),
  threshold_15_deaths = meta_long(rbind(cox_main, cox_extra[analysis == "extra_threshold_15"]))
)
meta_sens <- rbindlist(sens, idcol = "analysis")

# leave one dataset out
loo <- rbindlist(lapply(sets_main, function(nm) {
  m <- meta_long(cox_main[dataset != nm])
  m[, left_out_dataset := nm]
  m
}))

# ============================================================
# 8) Group tests
# ============================================================
# Proliferation genes: Hallmark E2F targets and G2M checkpoint (MSigDB)
h <- as.data.table(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"))
prolif <- unique(h[gs_name %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"), gene_symbol])
# MSigDB uses current symbols and the arrays those of their time: the proliferation
# pool is the part of the list that each platform recognizes
msg("[INFO] proliferation genes (E2F + G2M): %d | candidates among them: %d of %d (%s)",
    length(prolif), sum(genes %in% prolif), length(genes),
    paste(intersect(genes, prolif), collapse = ", "))

# Three nulls for the same question, on the mean of the meta-analytic logHR
# (random effects) of the candidates and on the fraction with HR < 1:
#   a) random gene sets, b) random proliferation genes,
#   c) permutation of survival within each dataset.
group_test <- function(data, cox_cand, label) {
  msg("[INFO] (%s) genome-wide Cox in %d datasets...", label, length(data))
  gw <- rbindlist(lapply(data, cox_genome))

  # the fast fit must match coxph on the candidates
  chk <- merge(cox_cand[, .(dataset, gene, b1 = logHR, s1 = se)],
               gw[, .(dataset, gene, b2 = logHR, s2 = se)], by = c("dataset", "gene"))
  stopifnot(nrow(chk) == nrow(cox_cand),
            max(abs(chk$b1 - chk$b2)) < 1e-6, max(abs(chk$s1 - chk$s2)) < 1e-6)
  msg("[OK] (%s) fast fit = coxph in the %d models of the candidates", label, nrow(chk))

  Bg <- dcast(gw, gene ~ dataset, value.var = "logHR")
  Sg <- dcast(gw, gene ~ dataset, value.var = "se")
  gene_univ <- Bg$gene
  Bg <- as.matrix(Bg[, -1]); Sg <- as.matrix(Sg[, -1])
  rownames(Bg) <- rownames(Sg) <- gene_univ
  Sg[is.na(Bg)] <- NA

  meta_c <- meta_long(cox_cand)
  g_eval <- meta_c$gene
  pattern_id <- vapply(g_eval, function(g) paste(colnames(Bg)[!is.na(Bg[g, ])], collapse = "|"),
                       character(1))

  # Meta-analysis of every gene measured in the same datasets as each candidate,
  # restricted to those datasets: each drawn gene has the same platform
  # coverage as the candidate it replaces
  pool <- lapply(unique(pattern_id), function(pid) {
    cols <- strsplit(pid, "|", fixed = TRUE)[[1]]
    ok   <- rowSums(is.na(Bg[, cols, drop = FALSE])) == 0
    m    <- meta_mat(Bg[ok, cols, drop = FALSE], Sg[ok, cols, drop = FALSE])
    data.table(gene = gene_univ[ok], logHR_random = m$logHR_random)
  })
  names(pool) <- unique(pattern_id)

  obs_vec <- setNames(meta_c$logHR_random, meta_c$gene)[g_eval]
  obs <- c(mean_logHR = mean(obs_vec), frac_HR_below_1 = mean(obs_vec < 0))
  for (g in g_eval) stopifnot(abs(pool[[pattern_id[g]]][gene == g, logHR_random] - obs_vec[g]) < 1e-8)

  # draw of distinct genes, one per candidate and with its coverage pattern
  draw <- function(allowed, n_rep) {
    idx <- lapply(g_eval, function(g) {
      p <- pool[[pattern_id[g]]]
      p[gene %in% allowed & !gene %in% genes]
    })
    size <- vapply(idx, nrow, integer(1))
    res <- mclapply(seq_len(n_rep), function(r) {
      repeat {
        pick <- vapply(seq_along(idx), function(i) sample.int(size[i], 1), integer(1))
        gs   <- vapply(seq_along(idx), function(i) idx[[i]]$gene[pick[i]], character(1))
        if (!anyDuplicated(gs)) break
      }
      v <- vapply(seq_along(idx), function(i) idx[[i]]$logHR_random[pick[i]], numeric(1))
      c(mean(v), mean(v < 0))
    }, mc.cores = N_CORES, mc.set.seed = TRUE)
    list(null = do.call(rbind, res), size = setNames(size, g_eval))
  }

  summarize_null <- function(null, name, pool_size) {
    rbindlist(lapply(1:2, function(j) {
      x <- null[, j]; o <- obs[j]
      data.table(
        analysis = label, test = name, statistic = names(obs)[j], observed = o,
        null_mean = mean(x), null_sd = sd(x),
        null_p2.5 = quantile(x, 0.025), null_p97.5 = quantile(x, 0.975),
        p_two_sided = (1 + sum(abs(x - mean(x)) >= abs(o - mean(x)))) / (length(x) + 1),
        p_protective_tail = if (j == 1) (1 + sum(x <= o)) / (length(x) + 1)
                            else (1 + sum(x >= o)) / (length(x) + 1),
        n_replicates = length(x), pool_size_min = min(pool_size), pool_size_max = max(pool_size))
    }))
  }

  set.seed(SEED)
  rnd <- draw(gene_univ, N_RAND)
  set.seed(SEED)
  pr <- draw(intersect(gene_univ, prolif), N_RAND)
  msg("[INFO] (%s) pool size by coverage pattern, random genes: %s", label,
      paste(sprintf("%s=%d", names(rnd$size), rnd$size)[!duplicated(pattern_id)], collapse = ", "))
  msg("[INFO] (%s) pool size by coverage pattern, proliferation genes: %s", label,
      paste(sprintf("%s=%d", names(pr$size), pr$size)[!duplicated(pattern_id)], collapse = ", "))

  # c) permuting survival keeps the correlation between the candidates,
  # which the gene draws do not respect
  Zs <- lapply(data, function(d) t(d$X[cox_cand[dataset == d$dataset, gene], , drop = FALSE]))
  perm_stat <- function(permute) {
    r <- rbindlist(lapply(names(data), function(nm) {
      d <- data[[nm]]
      i <- if (permute) sample.int(length(d$time)) else seq_along(d$time)
      f <- vapply(colnames(Zs[[nm]]),
                  function(g) cox_fast(Zs[[nm]][, g], d$time[i], d$event[i])[1:2], numeric(2))
      data.table(dataset = nm, gene = colnames(f), logHR = f[1, ], se = f[2, ])
    }))
    B <- dcast(r, gene ~ dataset, value.var = "logHR")
    S <- dcast(r, gene ~ dataset, value.var = "se")
    m <- meta_mat(as.matrix(B[, -1]), as.matrix(S[, -1]))$logHR_random
    c(mean(m), mean(m < 0))
  }
  stopifnot(max(abs(perm_stat(FALSE) - obs)) < 1e-8)
  set.seed(SEED)
  null_perm <- do.call(rbind, mclapply(seq_len(N_PERM), function(i) perm_stat(TRUE),
                                       mc.cores = N_CORES, mc.set.seed = TRUE))

  # where the set of all genes measured in every dataset falls
  genome <- meta_mat(Bg, Sg)
  genome[, gene := gene_univ]
  ref <- genome[k == length(data) & !gene %in% genes]
  msg("[INFO] (%s) genes measured in the %d datasets: %d | median logHR: %.4f | %% with HR < 1: %.1f",
      label, length(data), nrow(ref), median(ref$logHR_random), 100 * mean(ref$logHR_random < 0))

  list(group = rbind(summarize_null(rnd$null, "random_genes", rnd$size),
                     summarize_null(pr$null, "proliferation_genes_E2F_G2M", pr$size),
                     summarize_null(null_perm, "survival_permutation", NA_integer_)),
       nulls = list(random = rnd$null, prolif = pr$null, perm = null_perm), obs = obs, genome = genome)
}

pg <- list(main = group_test(data_main, cox_main, "main"))
for (sc in c("high_grade_2_3", "high_grade_package")) {
  cs <- cox_sub[analysis == sc]
  ds <- lapply(setNames(unique(cs$dataset), unique(cs$dataset)), prepare, subset = sc)
  pg[[sc]] <- group_test(ds, cs, sc)
}
group <- rbindlist(lapply(pg, `[[`, "group"))

# ============================================================
# 9) Comparison with TCGA
# ============================================================
tcga <- fread(tcga_tsv)
tcga <- tcga[, .(gene, HR_TCGA = HR_per_SD, HR_TCGA_lo95 = HR_lo95, HR_TCGA_hi95 = HR_hi95,
                 p_TCGA = p_value, BH_TCGA = p_bh)]
tcga[, `:=`(logHR_TCGA = log(HR_TCGA),
            se_TCGA = (log(HR_TCGA_hi95) - log(HR_TCGA_lo95)) / (2 * qnorm(0.975)))]

comp <- merge(tcga, meta_main[, .(gene, k, n_patients, n_deaths, HR_indep = HR_random,
                                  HR_indep_lo95 = HR_random_lo95, HR_indep_hi95 = HR_random_hi95,
                                  p_indep = p_random, BH_indep = BH_random, I2,
                                  logHR_indep = logHR_random, se_indep = se_random)],
              by = "gene", all = TRUE)
comp[, `:=`(
  same_direction = sign(logHR_TCGA) == sign(logHR_indep),
  # do the two estimates differ more than expected by chance?
  p_difference = 2 * pnorm(-abs((logHR_TCGA - logHR_indep) / sqrt(se_TCGA^2 + se_indep^2))))]
# the high-grade sensitivity analysis alongside: TCGA is almost entirely grade 2-3
hg <- sens$high_grade_2_3[, .(gene, k_HG23 = k, HR_indep_HG23 = HR_random, HR_indep_HG23_lo95 = HR_random_lo95,
                              HR_indep_HG23_hi95 = HR_random_hi95, p_indep_HG23 = p_random,
                              BH_indep_HG23 = BH_random, I2_HG23 = I2)]
comp <- merge(comp, hg, by = "gene", all.x = TRUE)
setorder(comp, p_TCGA)

rho <- cor.test(comp$logHR_TCGA, comp$logHR_indep, method = "spearman", exact = FALSE)

# ============================================================
# 10) Outputs
# ============================================================
meta_out <- merge(meta_main, tcga[, .(gene, HR_TCGA, HR_TCGA_lo95, HR_TCGA_hi95, p_TCGA, BH_TCGA)],
                  by = "gene", all = TRUE)
setorder(meta_out, p_random)
setcolorder(meta_out, c("gene", "k", "n_patients", "n_deaths", "n_HR_below_1",
                        "HR_random", "HR_random_lo95", "HR_random_hi95", "p_random", "BH_random",
                        "HR_fixed", "HR_fixed_lo95", "HR_fixed_hi95", "p_fixed", "BH_fixed",
                        "I2", "tau2", "Q", "p_Q"))

alias_tab <- unique(cox_all[, .(gene, dataset, feature_used)])
alias_tab <- dcast(alias_tab, gene ~ sub("_eset$", "", dataset), value.var = "feature_used", fill = "absent")

fwrite(dsets,     file.path(out_dir, "datasets_used.tsv"), sep = "\t")
fwrite(cox_all,   file.path(out_dir, "cox_by_dataset.tsv"), sep = "\t")
fwrite(meta_out,  file.path(out_dir, "meta_by_gene.tsv"), sep = "\t")
fwrite(meta_sens, file.path(out_dir, "meta_sensitivity.tsv"), sep = "\t")
fwrite(loo,       file.path(out_dir, "leave_one_out.tsv"), sep = "\t")
fwrite(group,     file.path(out_dir, "group_tests.tsv"), sep = "\t")
fwrite(comp,      file.path(out_dir, "comparison_TCGA.tsv"), sep = "\t")
fwrite(alias_tab, file.path(out_dir, "aliases_used.tsv"), sep = "\t")
saveRDS(list(datasets = dsets, cox = cox_all, meta = meta_out, sens = meta_sens, loo = loo,
             group = group, comp = comp, nulls = lapply(pg, `[[`, "nulls"),
             obs = lapply(pg, `[[`, "obs"), genome_meta = pg$main$genome, rho = rho),
        file.path(out_dir, "results.rds"))

# ============================================================
# 11) Summary in the log
# ============================================================
r3 <- function(d) { d <- copy(d); num <- names(d)[vapply(d, is.numeric, logical(1))]
  d[, (num) := lapply(.SD, signif, 3), .SDcols = num]; as.data.frame(d) }

cat("\n===== DATASETS =====\n")
print(as.data.frame(dsets[, .(dataset, platform, n_blacklisted, n_serous_OS, deaths_serous,
                              candidates_measured, status, reason)]))
cat("\n[TOTAL main] patients:", dsets[status == "main", sum(n_serous_OS)],
    "| deaths:", dsets[status == "main", sum(deaths_serous)], "\n")

cat("\n===== META-ANALYSIS PER GENE (main) =====\n")
print(r3(meta_out[, .(gene, k, n_patients, n_deaths, n_HR_below_1, HR_random, HR_random_lo95,
                      HR_random_hi95, p_random, BH_random, HR_fixed, p_fixed, BH_fixed, I2, HR_TCGA, p_TCGA)]))
cat("\n[RESULT] nominal p < 0.05 (random effects):", sum(meta_out$p_random < 0.05, na.rm = TRUE),
    "| BH < 0.05:", sum(meta_out$BH_random < 0.05, na.rm = TRUE),
    "| fixed effects nominal:", sum(meta_out$p_fixed < 0.05, na.rm = TRUE),
    "| fixed effects BH < 0.05:", sum(meta_out$BH_fixed < 0.05, na.rm = TRUE),
    "| HR < 1:", sum(meta_out$HR_random < 1, na.rm = TRUE), "of", sum(!is.na(meta_out$HR_random)), "\n")

cat("\n===== SENSITIVITY ANALYSES (summary) =====\n")
print(as.data.frame(meta_sens[, .(genes = .N, k_max = max(k), n_max = max(n_patients),
                                  deaths_max = max(n_deaths),
                                  nominal = sum(p_random < 0.05), BH_sig = sum(BH_random < 0.05),
                                  HR_below_1 = sum(HR_random < 1),
                                  mean_logHR = signif(mean(logHR_random), 3),
                                  min_p = signif(min(p_random), 3), gene_min_p = gene[which.min(p_random)]),
                              by = analysis]))

# Detail of the top gene of the discovery cohort and of the top gene of the meta-analysis
top_tcga <- tcga[which.min(p_TCGA), gene]
top_meta <- meta_main[which.min(p_random), gene]

cat(sprintf("\n===== %s (lowest p in TCGA) =====\n", top_tcga))
print(r3(cox_all[gene == top_tcga & analysis == "main",
                 .(dataset, n, deaths, HR, HR_lo95, HR_hi95, p, adjustment_covariates, HR_adj, p_adj)]))
print(r3(meta_sens[gene == top_tcga, .(analysis, k, n_patients, n_deaths, HR_random, HR_random_lo95,
                                       HR_random_hi95, p_random, HR_fixed, p_fixed, I2)]))
print(r3(loo[gene == top_tcga, .(left_out_dataset, k, HR_random, HR_random_lo95, HR_random_hi95, p_random)]))
gm <- meta_main[gene == top_tcga]
if (nrow(gm))
  cat(sprintf("[%s] one-sided p in the direction of TCGA (HR %s 1): %.3f\n", top_tcga,
              ifelse(tcga[gene == top_tcga, logHR_TCGA] < 0, "<", ">"),
              pnorm(-sign(tcga[gene == top_tcga, logHR_TCGA]) * gm$logHR_random / gm$se_random)))

cat(sprintf("\n===== %s (lowest p in the meta-analysis), leaving one dataset out =====\n", top_meta))
print(r3(loo[gene == top_meta, .(left_out_dataset, k, HR_random, HR_random_lo95, HR_random_hi95, p_random, BH_random)]))

cat("\n===== GROUP TESTS =====\n")
print(r3(group))
cat("\n[LOO] mean logHR of the candidates when each dataset is left out:\n")
print(as.data.frame(loo[, .(mean_logHR = signif(mean(logHR_random), 3), HR_below_1 = sum(logHR_random < 0),
                            nominal = sum(p_random < 0.05)), by = left_out_dataset]))

cat("\n===== COMPARISON WITH TCGA =====\n")
print(r3(comp[, .(gene, HR_TCGA, p_TCGA, HR_indep, HR_indep_lo95, HR_indep_hi95, p_indep, BH_indep,
                  k, same_direction, p_difference, HR_indep_HG23, p_indep_HG23, BH_indep_HG23)]))
rho_hg <- cor.test(comp$logHR_TCGA, log(comp$HR_indep_HG23), method = "spearman", exact = FALSE)
cat(sprintf("[RESULT] high grade 2-3: same direction as TCGA: %d of %d | Spearman: rho = %.3f, p = %.3g\n",
            sum(sign(comp$logHR_TCGA) == sign(log(comp$HR_indep_HG23))), nrow(comp),
            rho_hg$estimate, rho_hg$p.value))
cat(sprintf("\n[RESULT] same direction as TCGA: %d of %d | Spearman logHR TCGA vs independent: rho = %.3f, p = %.3g\n",
            sum(comp$same_direction, na.rm = TRUE), sum(!is.na(comp$same_direction)),
            rho$estimate, rho$p.value))

cat("\n===== SESSION =====\n")
print(sessionInfo())
