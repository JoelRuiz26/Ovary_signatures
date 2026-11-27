#!/usr/bin/env Rscript

# ============================================================
# Fallopian tube contrast, step 1: download the candidate GEO series
# (GPL570, tumor and fallopian tube epithelium) and dump their phenotype
# metadata, without assuming anything about the groups.
#
# Output (9_Revision_analyses/fallopian_tube/):
#         geo_data/<GSE>_series_matrix.txt.gz   downloaded with GEOquery
#         pheno_<GSE>_<GPL>.tsv, pheno_candidate_series.rds
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})
options(timeout = 600)

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/fallopian_tube")
GEO_DIR <- file.path(OUT_DIR, "geo_data")
dir.create(GEO_DIR, recursive = TRUE, showWarnings = FALSE)

series <- c("GSE10971", "GSE69428", "GSE107927", "GSE52037", "GSE28044")
res <- list()
for (g in series) {
  cat("\n\n==================== ", g, " ====================\n")
  es <- tryCatch(getGEO(g, GSEMatrix = TRUE, getGPL = FALSE, destdir = GEO_DIR),
                 error = function(e) {cat("ERROR:", conditionMessage(e), "\n"); NULL})
  if (is.null(es)) next
  for (i in seq_along(es)) {
    e  <- es[[i]]
    pd <- pData(e)
    cat("-- ExpressionSet", i, "| platform:", annotation(e), "| samples:", ncol(e),
        "| probes:", nrow(e), "\n")
    cat("Series title:", experimentData(e)@title, "\n")
    rng <- range(exprs(e), na.rm = TRUE)
    cat("Expression range:", signif(rng, 4), "\n")
    cat("Processing:", unique(pd$data_processing)[1:min(3, length(unique(pd$data_processing)))], "\n")
    # informative columns
    cols <- colnames(pd)[vapply(pd, function(x) {n <- length(unique(as.character(x))); n > 1 && n <= 15}, logical(1))]
    for (cc in cols) {cat("\n---", cc, "---\n"); print(table(pd[[cc]], useNA = "ifany"))}
    keep <- intersect(c("title", "geo_accession", "source_name_ch1",
                        grep(":ch1$|^characteristics_ch1", colnames(pd), value = TRUE), "description"), colnames(pd))
    out <- pd[, keep, drop = FALSE]
    write.table(out, file.path(OUT_DIR, sprintf("pheno_%s_%s.tsv", g, annotation(e))),
                sep = "\t", quote = FALSE, row.names = FALSE)
    res[[paste(g, annotation(e))]] <- out
  }
}
saveRDS(res, file.path(OUT_DIR, "pheno_candidate_series.rds"))
sessionInfo()
