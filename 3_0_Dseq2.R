suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
})

## =========================
## 0) Paths & settings
## =========================
DDSDIR    <- "~/Ovary_signatures/2_1_Normalized_corrected_counts/"
OUTDIR_DE <- "~/Ovary_signatures/3_1_DE_results/"
dir.create(OUTDIR_DE, showWarnings = FALSE, recursive = TRUE)

LOG2FC_TH <- 1
PADJ_TH   <- 0.01

dds_info <- list(
  auto = list(
    file   = "dds_tumor_CtlAuto.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_autoencoder"
  ),
  homolog = list(
    file   = "dds_tumor_Ctlhomolog.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_GTEx"
  ),
  clean_no_neg = list(
    file   = "dds_tumor_Ctlhomolog_CLEAN_NON_NEGATIVE.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_GTEx"
  ),
  clean_rescued = list(
    file   = "dds_tumor_Ctlhomolog_CLEAN_rescuedNegatives.rds",
    factor = "condition",
    num    = "TCGA_tumor",
    den    = "control_GTEx"
  )
)

## =========================
## 1) Helper: get results safely
## =========================
get_results_safe <- function(dds, contrast_vec) {
  # Primero intentamos results()
  res_try <- try(results(dds, contrast = contrast_vec), silent = TRUE)
  
  if (!inherits(res_try, "try-error")) {
    return(list(dds = dds, res = res_try, refit = FALSE))
  }
  
  # Si falla, corremos DESeq (una vez) y reintentamos
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = contrast_vec)
  
  list(dds = dds, res = res, refit = TRUE)
}

## =========================
## 2) Run DE
## =========================
run_DE <- function(tag, info) {
  message("\n========== ", tag, " ==========")
  
  dds_path <- file.path(DDSDIR, info$file)
  if (!file.exists(dds_path)) stop("No existe: ", dds_path)
  
  dds <- readRDS(dds_path)
  
  # Chequeo mínimo de factor/levels
  if (!(info$factor %in% colnames(colData(dds)))) {
    stop("No encuentro '", info$factor, "' en colData para ", info$file)
  }
  lev <- levels(colData(dds)[[info$factor]])
  if (!all(c(info$num, info$den) %in% lev)) {
    stop("Niveles faltantes en ", info$file, ": ",
         info$num, " o ", info$den, " (levels=", paste(lev, collapse=", "), ")")
  }
  
  contrast_vec <- c(info$factor, info$num, info$den)
  out <- get_results_safe(dds, contrast_vec)
  dds_used <- out$dds
  res <- out$res
  
  if (out$refit) {
    message("-> Se corrió DESeq() porque el objeto no estaba completamente ajustado.")
  } else {
    message("-> Se usó el ajuste existente del objeto (no se re-corrió DESeq()).")
  }
  
  res <- res[order(res$pvalue), ]
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df %>% dplyr::select(gene, dplyr::everything())
  
  # Conteos de Up/Down/NS para sanity-check rápido
  res_df <- res_df %>%
    mutate(
      Type = case_when(
        !is.na(padj) & padj < PADJ_TH & log2FoldChange >=  LOG2FC_TH ~ "Upregulated",
        !is.na(padj) & padj < PADJ_TH & log2FoldChange <= -LOG2FC_TH ~ "Downregulated",
        TRUE ~ "NotSignificant"
      )
    )
  
  tab <- table(res_df$Type)
  message("Up: ", tab[["Upregulated"]] %||% 0,
          " | Down: ", tab[["Downregulated"]] %||% 0,
          " | NS: ", tab[["NotSignificant"]] %||% 0,
          " | Total: ", nrow(res_df))
  
  res_DEG <- res_df %>%
    filter(Type %in% c("Upregulated", "Downregulated"))
  
  # Save FULL + DEG
  write.table(res_df,
              file = file.path(OUTDIR_DE, paste0(tag, "_FULL.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(res_DEG,
              file = file.path(OUTDIR_DE, paste0(tag, "_DEG.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Guardar el dds realmente usado (ya ajustado si fue necesario)
  saveRDS(dds_used,
          file.path(OUTDIR_DE, paste0(tag, "_dds_used.rds")))
  
  invisible(list(dds = dds_used, full = res_df, deg = res_DEG))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

## =========================
## 3) Run all 4
## =========================
out_list <- list()
for (nm in names(dds_info)) {
  out_list[[nm]] <- run_DE(nm, dds_info[[nm]])
}
