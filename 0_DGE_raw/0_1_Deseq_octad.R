#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
})

# ===================== CONFIG ===================== #
out_dir <- "~/Ovary_signatures/0_DGE_raw/"
setwd(out_dir)

h5_path <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"

lfc_thr  <- 1.0
padj_thr <- 0.01

# ===================== METADATA ===================== #
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

# Cases: TCGA OVARY (no adjacent)
case_ovary <- phenoDF %>%
  filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
         data.source == "TCGA",
         sample.type != "adjacent") %>%
  pull(sample.id)

# Controls: GTEx normal + TCGA adjacent OVARY
controls_all <- phenoDF %>%
  filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
         (data.source == "GTEX" & sample.type == "normal") |
           (data.source == "TCGA" & sample.type == "adjacent")) %>%
  pull(sample.id)

metadata <- phenoDF %>% filter(sample.id %in% c(case_ovary, controls_all))
saveRDS(metadata, "0_1_metadata.rds")

cat("Cases:", length(case_ovary), "| Controls:", length(controls_all), "\n")

# ===================== RUN DESeq2 via OCTAD ===================== #
res <- octad::diffExp(
  case_id           = case_ovary,
  control_id        = controls_all,
  source            = "octad.whole",
  file              = h5_path,
  normalize_samples = TRUE,   # OCTAD usa RUVg para DESeq2 cuando TRUE
  k                 = 2,      # como ya lo traías
  n_topGenes        = 15000,
  DE_method         = "DESeq2",
  annotate          = TRUE,
  output            = FALSE
)

# ---- FIX: evita pvalue/padj = 0 por underflow ----
eps <- .Machine$double.xmin
res$pvalue_safe <- if ("pvalue" %in% names(res)) pmax(res$pvalue, eps) else NA_real_
res$padj_safe   <- if ("padj"   %in% names(res)) pmax(res$padj,   eps) else NA_real_
res$mlog10_p    <- -log10(res$pvalue_safe)
res$mlog10_padj <- -log10(res$padj_safe)

# (opcional) reemplazar en sitio para que NUNCA queden 0 en columnas originales
res$pvalue <- res$pvalue_safe
res$padj   <- res$padj_safe

# Guardar resultado completo
saveRDS(res, "DE_full_OVARY_DESeq2.rds")

# ===================== BUILD SIGNATURES (ALL + SIGNIFIC) ===================== #
# Usa Symbol si existe; si no, cae a identifier
if ("Symbol" %in% names(res)) {
  res2 <- res %>% mutate(Symbol = toupper(Symbol))
} else {
  res2 <- res %>% mutate(Symbol = identifier)
}

sig_all <- res2 %>%
  filter(!is.na(log2FoldChange), !is.na(pvalue), !is.na(padj)) %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  transmute(
    Symbol,
    log2FoldChange,
    regulation = if_else(log2FoldChange > 0, "up", "down"),
    pvalue,
    padj
  )

sig_sig <- sig_all %>%
  filter(abs(log2FoldChange) >= lfc_thr, padj <= padj_thr)

saveRDS(sig_all, "DE_OVARY_ALL_DESeq2.rds")
saveRDS(sig_sig, "DE_OVARY_signific_DESeq2.rds")

# ===================== QUICK QC PRINT ===================== #
cnt <- function(x) {
  tb <- table(x)
  up   <- if ("up"   %in% names(tb)) tb[["up"]]   else 0L
  down <- if ("down" %in% names(tb)) tb[["down"]] else 0L
  c(up = as.integer(up), down = as.integer(down))
}

c_all <- cnt(sig_all$regulation)
c_sig <- cnt(sig_sig$regulation)

cat("OVARY | DESeq2",
    "-> TOTAL:", nrow(sig_all), "(UP:", c_all["up"], "| DOWN:", c_all["down"], ")",
    "| SIGNIFIC:", nrow(sig_sig), "(UP:", c_sig["up"], "| DOWN:", c_sig["down"], ")\n")

cat("Zeros check (should be 0): pvalue==0:", sum(res$pvalue == 0, na.rm = TRUE),
    "| padj==0:", sum(res$padj == 0, na.rm = TRUE), "\n")

save.image("Image_DiseaseSignature_OVARY_DESeq2_only.RData")
cat("Done.\n")
