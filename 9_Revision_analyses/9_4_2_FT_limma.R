#!/usr/bin/env Rscript

# ============================================================
# Fallopian tube contrast, step 2: limma of tumor (HGSC) against normal
# fallopian tube epithelium, in each GPL570 series.
#
# Same procedure as 2_DEG_GEO/GEO_multi_RRA.R: series matrix, log2(x + 1)
# if the data are linear, lmFit + eBayes, frozen probe -> symbol
# dictionary and one probe per gene (largest B). The two series are
# combined with Stouffer's method weighted by the square root of the
# number of arrays.
#
# Input:  9_Revision_analyses/fallopian_tube/geo_data/GSE10971_series_matrix.txt.gz,
#         .../GSE69428_series_matrix.txt.gz (9_4_1_FT_series_metadata.R),
#         2_DEG_GEO/annotation/gprofiler_GPL570_probe_to_hgnc.tsv.gz (repository)
# Output (9_Revision_analyses/fallopian_tube/):
#         de_<contrast>.tsv, de_Stouffer.tsv, de_all.rds,
#         samples_by_series.tsv, series_used.tsv
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(dplyr)
})
set.seed(1)

# Folder of this script, to find the frozen probe annotation of the repository
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir  <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
PROBE_DICT  <- Sys.getenv("GEO_PROBE_DICT",
                          file.path(script_dir, "../2_DEG_GEO/annotation/gprofiler_GPL570_probe_to_hgnc.tsv.gz"))

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/fallopian_tube")
GEO_DIR <- file.path(OUT_DIR, "geo_data")

if (!file.exists(PROBE_DICT))
  stop("Probe dictionary not found: ", PROBE_DICT, " (2_DEG_GEO/GEO_multi_RRA.R builds it; or set GEO_PROBE_DICT)")
dict <- read.delim(PROBE_DICT, comment.char = "#", colClasses = "character", na.strings = "NA")

# --- helpers (same logic as the main pipeline) ---
log2_auto <- function(expr) {
  rng <- range(expr, na.rm = TRUE)
  if (rng[2] > 100 || (rng[2] - rng[1]) > 50) {
    message("   linear scale -> log2(expr + 1)"); log2(expr + 1)
  } else {message("   already in log2"); expr}
}
collapse_probes <- function(tt) {
  stopifnot(all(rownames(tt) %in% dict$probe))
  tt$Gene.symbol <- dict$hgnc_symbol[match(rownames(tt), dict$probe)]
  tt <- tt[!is.na(tt$Gene.symbol) & tt$Gene.symbol != "", ]
  tt <- tt[order(tt$Gene.symbol, -tt$B), ]
  tt <- tt[!duplicated(tt$Gene.symbol), ]
  rownames(tt) <- tt$Gene.symbol
  tt
}
# limma fit; the design must have the coefficient "groupcase"
fit_contrast <- function(expr, design, name) {
  fit <- eBayes(lmFit(expr, design))
  tt <- topTable(fit, coef = "groupcase", number = Inf, adjust.method = "BH", sort.by = "P")
  tt$probe <- rownames(tt)
  tt$df_total <- setNames(fit$df.total, rownames(fit))[rownames(tt)]
  tt <- collapse_probes(tt)
  # adjusted p recomputed at the gene level (one probe per gene)
  tt$adj.P.gene <- p.adjust(tt$P.Value, "BH")
  tt$z <- zscoreT(tt$t, tt$df_total)
  stopifnot(all(is.finite(tt$z)))
  out <- tt[, c("Gene.symbol", "probe", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "adj.P.gene", "B", "z")]
  colnames(out)[1] <- "gene"
  message(sprintf("   %s: %d genes | adj.P<0.05 (probe level, as in the pipeline): %d up, %d down",
                  name, nrow(out), sum(out$adj.P.Val < 0.05 & out$logFC > 0), sum(out$adj.P.Val < 0.05 & out$logFC < 0)))
  write.table(out, file.path(OUT_DIR, sprintf("de_%s.tsv", name)), sep = "\t", quote = FALSE, row.names = FALSE)
  out
}
# secretory / ciliated epithelium, proliferation and stroma markers, printed as a sanity check
markers <- c("PAX8", "WT1", "MKI67", "TOP2A", "CDKN2A", "FOXJ1", "CAPS", "OVGP1", "EPCAM", "KRT8",
             "DCN", "COL1A1", "ACTA2", "PTPRC", "ZEB1", "PRDM6", "OSR1", "LHX9")
show_markers <- function(de) print(round(de[intersect(markers, rownames(de)), c("logFC", "t", "adj.P.Val")], 3))

samples <- list(); de <- list(); n_arr <- c()

# =====================================================================
# GSE10971 (Tone et al. 2008): laser-capture microdissected distal FTE and tubal/ovarian HGSC
# =====================================================================
message("\n=== GSE10971 ===")
e1 <- getGEO(filename = file.path(GEO_DIR, "GSE10971_series_matrix.txt.gz"), getGPL = FALSE)
p1 <- pData(e1)
p1$group0 <- dplyr::recode(p1$source_name_ch1,
  "Laser-capture microdissected high grade serous carcinoma" = "HGSC",
  "Laser-capture microdissected non-malignant FTE from normal control patient" = "FTEn",
  "Laser-capture microdissected non-malignant FTE from BRCA1/2-mutation carrier" = "FTEb")
p1$site <- ifelse(p1$group0 == "HGSC", ifelse(grepl("of ovary", p1$description), "ovary", "tube"), NA)
# patient: among the tumors, patients 2 and 3 have two arrays (tubal and ovarian tumor)
p1$patient <- paste0(p1$group0, "_", sub(".*patient ([0-9]+).*", "\\1", p1$title))
print(table(p1$group0, p1$site, useNA = "ifany"))
print(p1[p1$group0 == "HGSC", c("title", "site", "patient", "BRCA1/2 mutation:ch1")])
dup <- p1$patient[duplicated(p1$patient)]
message("Patients with two tumor arrays: ", paste(dup, collapse = ", "))
x1 <- log2_auto(exprs(e1))

# Main contrast: one array per patient (for patients 2 and 3 the tumor of ovarian site is kept)
# against FTE of controls without mutation (FTEn). FTEb (BRCA1/2 carriers) is analyzed separately.
excl_dup <- p1$patient %in% dup & p1$site == "tube"
k <- (p1$group0 == "HGSC" & !excl_dup) | p1$group0 == "FTEn"
g <- factor(ifelse(p1$group0[k] == "HGSC", "case", "control"), levels = c("control", "case"))
print(table(g))
n_arr["GSE10971"] <- sum(k)
de$GSE10971 <- fit_contrast(x1[, k], model.matrix(~ group, data.frame(group = g)), "GSE10971")
show_markers(de$GSE10971)

# Sensitivity a: the 13 tumor arrays (as the authors) against FTEn
k <- p1$group0 %in% c("HGSC", "FTEn")
g <- factor(ifelse(p1$group0[k] == "HGSC", "case", "control"), levels = c("control", "case"))
de$GSE10971_13tum <- fit_contrast(x1[, k], model.matrix(~ group, data.frame(group = g)), "GSE10971_13tum")
# Sensitivity b: 11 tumors against FTEn + FTEb (24 epithelia)
k <- (p1$group0 == "HGSC" & !excl_dup) | p1$group0 %in% c("FTEn", "FTEb")
g <- factor(ifelse(p1$group0[k] == "HGSC", "case", "control"), levels = c("control", "case"))
de$GSE10971_FTEn_FTEb <- fit_contrast(x1[, k], model.matrix(~ group, data.frame(group = g)), "GSE10971_FTEn_FTEb")
# Sensitivity c: tumors of ovarian site only (7) against FTEn
k <- (p1$group0 == "HGSC" & p1$site == "ovary") | p1$group0 == "FTEn"
g <- factor(ifelse(p1$group0[k] == "HGSC", "case", "control"), levels = c("control", "case"))
de$GSE10971_ovarian_site <- fit_contrast(x1[, k], model.matrix(~ group, data.frame(group = g)), "GSE10971_ovarian_site")
# Diagnostic: FTEb against FTEn (how different the epithelium of carriers is)
k <- p1$group0 %in% c("FTEn", "FTEb")
g <- factor(ifelse(p1$group0[k] == "FTEb", "case", "control"), levels = c("control", "case"))
d_b <- fit_contrast(x1[, k], model.matrix(~ group, data.frame(group = g)), "GSE10971_diag_FTEb_vs_FTEn")
# Diagnostic: luteal against follicular phase within FTEn
k <- p1$group0 == "FTEn"
g <- factor(ifelse(p1$`Cycle:ch1`[k] == "Luteal", "case", "control"), levels = c("control", "case"))
d_c <- fit_contrast(x1[, k], model.matrix(~ group, data.frame(group = g)), "GSE10971_diag_luteal_vs_follicular")

samples$GSE10971 <- data.frame(
  series = "GSE10971", gsm = p1$geo_accession, title = p1$title, group = p1$group0, patient = p1$patient,
  detail = ifelse(p1$group0 == "HGSC", paste0("site ", p1$site, "; BRCA ", p1$`BRCA1/2 mutation:ch1`),
                  paste0("phase ", p1$`Cycle:ch1`, "; BRCA ", p1$`BRCA1/2 mutation:ch1`)),
  main_use = ifelse(p1$group0 == "FTEn", "control", ifelse(p1$group0 == "HGSC" & !excl_dup, "tumor", "excluded")),
  exclusion_reason = ifelse(p1$group0 == "FTEb", "non-malignant FTE from a BRCA1/2 mutation carrier (sensitivity analysis only)",
                     ifelse(excl_dup, "second array of the same patient (tubal tumor; the ovarian one is kept)", "")))

# =====================================================================
# GSE69428 (Yamamoto et al. 2016): HGSOC and paired normal oviduct, LCM, 10 patients
# =====================================================================
message("\n=== GSE69428 ===")
e2 <- getGEO(filename = file.path(GEO_DIR, "GSE69428_series_matrix.txt.gz"), getGPL = FALSE)
p2 <- pData(e2)
p2$group0 <- dplyr::recode(p2$`tissue:ch1`, "LCM from cryosection (HGSOC)" = "HGSOC",
                           "LCM from cryosection (normal oviduct)" = "normal_oviduct", .default = "other")
p2$patient <- ifelse(p2$group0 == "other", NA, sub("^(patient_[0-9]+)_.*", "\\1", p2$title))
print(table(p2$`tissue:ch1`, p2$group0))
print(table(p2$patient, p2$group0)[, c("HGSOC", "normal_oviduct")])
x2 <- log2_auto(exprs(e2))
k <- p2$group0 != "other"
g <- factor(ifelse(p2$group0[k] == "HGSOC", "case", "control"), levels = c("control", "case"))
pat <- factor(p2$patient[k])
# QC of the pairing: tumor-normal correlation within and between patients
cc <- cor(x2[, k][, g == "case"], x2[, k][, g == "control"], method = "spearman")
dimnames(cc) <- list(pat[g == "case"], pat[g == "control"])
cc <- cc[order(rownames(cc)), order(colnames(cc))]
message(sprintf("   pairing QC: tumor-normal rho, same patient = %.3f (mean); different patient = %.3f",
                mean(diag(cc)), mean(cc[row(cc) != col(cc)])))
# Main contrast: paired design (~ patient + group)
n_arr["GSE69428"] <- sum(k)
de$GSE69428 <- fit_contrast(x2[, k], model.matrix(~ pat + group, data.frame(pat = pat, group = g)), "GSE69428")
show_markers(de$GSE69428)
# Sensitivity: unpaired (~ group), identical to the main pipeline
de$GSE69428_unpaired <- fit_contrast(x2[, k], model.matrix(~ group, data.frame(group = g)), "GSE69428_unpaired")

samples$GSE69428 <- data.frame(
  series = "GSE69428", gsm = p2$geo_accession, title = p2$title, group = p2$group0, patient = p2$patient,
  detail = p2$`tissue:ch1`,
  main_use = ifelse(p2$group0 == "HGSOC", "tumor", ifelse(p2$group0 == "normal_oviduct", "control", "excluded")),
  exclusion_reason = ifelse(p2$group0 == "other", "cultured, immortalized or transformed FTSC, or xenograft (not primary tissue)", ""))

# =====================================================================
# Combined signature (Stouffer weighted by the square root of the number of arrays), common genes
# =====================================================================
common <- intersect(de$GSE10971$gene, de$GSE69428$gene)
z1 <- de$GSE10971[common, "z"]; z2 <- de$GSE69428[common, "z"]
w <- sqrt(n_arr[c("GSE10971", "GSE69428")])
comb <- data.frame(gene = common, z_GSE10971 = z1, z_GSE69428 = z2,
                   logFC_GSE10971 = de$GSE10971[common, "logFC"], logFC_GSE69428 = de$GSE69428[common, "logFC"],
                   z = (w[1] * z1 + w[2] * z2) / sqrt(sum(w^2)))
comb$P.Value <- 2 * pnorm(-abs(comb$z))
comb$adj.P.gene <- p.adjust(comb$P.Value, "BH")
comb$same_sign <- sign(comb$logFC_GSE10971) == sign(comb$logFC_GSE69428)
write.table(comb, file.path(OUT_DIR, "de_Stouffer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("\nStouffer: %d common genes | arrays: %s | Spearman rho of z between series = %.3f | same logFC sign = %.1f%% | FDR<0.05: %d up, %d down",
                nrow(comb), paste(names(n_arr), n_arr, collapse = ", "), cor(z1, z2, method = "spearman"),
                100 * mean(comb$same_sign),
                sum(comb$adj.P.gene < 0.05 & comb$z > 0), sum(comb$adj.P.gene < 0.05 & comb$z < 0)))

# --- documentation tables ---
smp <- do.call(rbind, samples)
write.table(smp, file.path(OUT_DIR, "samples_by_series.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(smp$series, smp$main_use))
n_use <- function(s, use) sum(smp$series == s & smp$main_use == use)

series <- data.frame(
  series = c("GSE10971", "GSE69428", "GSE107927", "GSE52037", "GSE28044"),
  platform = c("GPL570", "GPL570", "GPL570 (expression) + GPL11154 (ChIP-seq)", "GPL570", "GPL570"),
  used = c("yes", "yes", "no", "no", "no"),
  n_tumor = c(n_use("GSE10971", "tumor"), n_use("GSE69428", "tumor"), NA, NA, NA),
  n_fallopian_tube = c(n_use("GSE10971", "control"), n_use("GSE69428", "control"), NA, NA, NA),
  phenotype_column = c("source_name_ch1 (group); description (tumor site); title (patient); BRCA1/2 mutation:ch1",
                       "tissue:ch1 (group); title (patient, for the pairing)", NA, "source_name_ch1", "source_name_ch1"),
  limma_design = c("~ group (unpaired)", "~ patient + group (paired)", NA, NA, NA),
  deposited_preprocessing = c("RMA (ArrayAssist), linear scale -> log2(x+1); 54,675 probes",
                              "Expression Console/Partek, log2; deposited matrix filtered to 33,362 probes", NA, NA, NA),
  n_genes = c(nrow(de$GSE10971), nrow(de$GSE69428), NA, NA, NA),
  excluded_samples = c(
    "12 FTE from BRCA1/2 mutation carriers (FTEb; sensitivity analysis only) and 2 tubal tumor arrays of patients 2 and 3, who also have an ovarian tumor (one array per patient is kept)",
    "9 arrays that are not primary tissue: cultured FTSC (2), immortalized (2), transformed (2) and xenografts (3)",
    "whole series: only 4 serous carcinomas on GPL570, no normal fallopian tube with expression data",
    "whole series: the control is ovarian surface epithelium (10), not fallopian tube; 10 primary tumors",
    "whole series: non-malignant FTE only (12 BRCA1+ and 12 BRCA1-), no tumors; probe-centered data"),
  note = c("distal FTE by LCM, 6 in luteal and 6 in follicular phase; tumors: 6 of tubal and 7 of ovarian site (11 patients)",
           "LCM of cryosections; tumor and normal oviduct from the same 10 patients",
           NA, "checked: it is GPL570, but it has no fallopian tube",
           "same group as GSE10971; cannot be combined without confounding batch and group"))
write.table(series, file.path(OUT_DIR, "series_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
saveRDS(de, file.path(OUT_DIR, "de_all.rds"))
sessionInfo()
