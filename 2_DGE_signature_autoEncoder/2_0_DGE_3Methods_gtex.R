# =========================================
# Load Required Libraries
# =========================================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)
library(harmonicmeanp)  


setwd("~/Ovary_signatures/2_DGE_signature_autoEncoder/")

# =========================================
# Download HDF5
# =========================================
h5_path <- file.path("/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5")
#if (!file.exists(h5_path)) {
#  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5"
#  message(">> Descargando HDF5 (~3 GB) a: ", h5_path)
#  download.file(url, destfile = h5_path, mode = "wb", quiet = FALSE)
#}
#stopifnot(file.exists(h5_path))  # ensure

# =========================================
# Metadata (ExperimentHub: EH7274)
# =========================================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
dim(phenoDF)

# =========================================
# Load samples 
# =========================================
##Option 1: When you know the signature
muestras_ovary <- readRDS("0_OSC_tcga_ids.RDS")
length(muestras_ovary) #[1] 427
# Ensure IDs existence available in OCTAD
case_ovary <- intersect(muestras_ovary, phenoDF$sample.id)
if (length(case_ovary) < length(muestras_ovary)) warning("Missing some samples in OVARY")
#setdiff(muestras_ovary, case_ovary) #[1] "TCGA-25-1870-01"

##Option 2: When you don't know the signature
#case_ovary <- phenoDF %>% 
#  dplyr::filter(grepl("OVARY", biopsy.site, ignore.case = TRUE),
#    data.source == "TCGA",
#    sample.type != "adjacent") %>% dplyr::pull(sample.id)
#length(case_ovary)
#[1] 426

# =========================================
# Control samples (GTEx y TCGA-adjacent)
# =========================================
# load vector ordenado de controles seleccionado por el autoencoder
controls_all <- computeRefTissue(
  case_id      = case_ovary,
  adjacent     = FALSE,      #if TRUE = normal GTEX + adjacent TCGA & TARGET
  source       = "octad",    
  control_size = 88) %>% as.character()

length(controls_all) #[1] 88
metadata <- phenoDF %>% dplyr::filter(sample.id %in% c(controls_all,case_ovary))
metadata_control <- metadata %>% filter(sample.type %in% c("normal", "adjacent"))
table(metadata_control$biopsy.site, metadata_control$data.source)

saveRDS(metadata,
        "2_1_Output_gtex_rds/2_1_metadata.rds")
table(metadata$subtype)

# sample.id list by subtype
#ovca_subtypes_list <- metadata %>%
#  dplyr::filter(sample.id %in% case_ovary) %>%
#  split(.$subtype) %>%
#  lapply(function(x) x$sample.id)
#saveRDS(ovca_subtypes_list,"OVCA_TCGA_subtypes_ids_list.rds")

# =========================================
# DGE with 3 methods (solo OVARY)
# =========================================
lfc_thr  <- 1.0
padj_thr <- 0.01

pick_id_col <- function(df){
  if ("identifier" %in% names(df)) return("identifier")
  if ("ENSEMBL"   %in% names(df)) return("ENSEMBL")
  stop("No gene ID column found (identifier/ENSEMBL).")
}
pick_symbol_col <- function(df){
  if ("Symbol" %in% names(df)) return("Symbol")
  if ("Symbol_autho" %in% names(df)) return("Symbol_autho")
  NA_character_
}

make_signature_all <- function(res_df){
  sym_col <- pick_symbol_col(res_df)
  df <- if (!is.na(sym_col)) {
    res_df %>% dplyr::mutate(Symbol = toupper(.data[[sym_col]]))
  } else {
    res_df %>% dplyr::mutate(Symbol = .data[[pick_id_col(res_df)]])
  }
  
  df %>%
    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      !is.na(pvalue)          # (OPCIONAL pero recomendable)
    ) %>%
    dplyr::distinct(Symbol, .keep_all = TRUE) %>%
    dplyr::transmute(
      Symbol,
      log2FoldChange,
      regulation = dplyr::if_else(log2FoldChange > 0, "up", "down"),
      pvalue,               
      padj                 
    )
}


make_signature_sig <- function(res_df, lfc_thr = 1, padj_thr = 0.01){
  make_signature_all(res_df) %>%
    dplyr::filter(abs(log2FoldChange) >= lfc_thr, padj <= padj_thr)
}

run_and_save <- function(clade_label, case_ids, ctrl_ids, method = c("edgeR","DESeq2","limma")){
  method <- match.arg(method)
  method_tag <- switch(method, edgeR="EdgeR", DESeq2="DESeq2", limma="limma")
  normalize_flag <- method %in% c("edgeR","DESeq2")
  
  res <- diffExp(
    case_id           = case_ids,
    control_id        = ctrl_ids,
    source            = "octad.whole",
    file              = h5_path,
    normalize_samples = normalize_flag,
    k                 = if (normalize_flag) 2 else 0,
    n_topGenes        = 10000,
    DE_method         = method,
    annotate          = TRUE,
    output            = FALSE
  )
  
  saveRDS(res, file.path("2_1_Output_gtex_rds/",
                         paste0("DE_full_", clade_label, "_", method_tag, ".rds")))
  
  sig_all <- make_signature_all(res)
  sig_sig <- make_signature_sig(res, lfc_thr = lfc_thr, padj_thr = padj_thr)
  
  saveRDS(sig_all, file.path("2_1_Output_gtex_rds/",
                             paste0("DE_", clade_label, "_ALL_", method_tag, ".rds")))
  saveRDS(sig_sig, file.path("2_1_Output_gtex_rds/",
                             paste0("DE_", clade_label, "_signific_", method_tag, ".rds")))
  
  cnt <- function(x) {
    tb   <- table(x)
    up   <- if ("up"   %in% names(tb)) tb[["up"]]   else 0L
    down <- if ("down" %in% names(tb)) tb[["down"]] else 0L
    c(up=as.integer(up), down=as.integer(down))
  }
  c_all <- cnt(sig_all$regulation)
  c_sig <- cnt(sig_sig$regulation)
  
  cat(clade_label, "|", method_tag,
      "-> TOTAL:",   nrow(sig_all), "(UP:", c_all["up"],   "| DOWN:", c_all["down"],   ")",
      "| SIGNIFIC:", nrow(sig_sig), "(UP:", c_sig["up"],  "| DOWN:", c_sig["down"], ")\n")
}

# Ejecutar solo para OVARY
for (m in c("edgeR","DESeq2","limma")) run_and_save("OVARY", case_ovary, controls_all, m)
cat("Done. Signatures saved ALL and signific per method for OVARY\n")



save.image("/STORAGE/csbig/jruiz/Ovary_data/2_Output_AutoEncoder/gtex/Image_DiseaseSignature_OVARY.RData")
