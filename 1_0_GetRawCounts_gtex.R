###############################################
## Round raw counts to integers + save outputs
###############################################

suppressPackageStartupMessages({
  library(octad)
  library(octad.db)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(edgeR)
  library(ggfortify)
  library(Rtsne)
  library(uwot)
})

options(stringsAsFactors = FALSE)

## =========================
## 0) Paths and variables
## =========================
#Tumor tissue to explore
tumor <- as.character("OVARY")

OUT_DIR <-  "~/Ovary_signatures/1_1_Counts_Raw/"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

H5_PATH    <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"
#if (!file.exists(H5_PATH)) {
#  url <- "https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tp>
#  message(">> Descargando HDF5 (~3 GB) a: ", H5_PATH)
#  download.file(url, destfile = H5_PATH, mode = "wb", quiet = FALSE)
#}
#stopifnot(file.exists(H5_PATH))  # ensure

## =========================
## 1) Load OCTAD metadata
## =========================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

## =========================
## 2) Define OVARY tumor samples (TCGA)
## =========================
case_tumor <- phenoDF %>% 
  filter(
    grepl(tumor, biopsy.site, ignore.case = TRUE),
    data.source == "TCGA",
    sample.type != "adjacent"
  ) %>% 
  pull(sample.id)

print(length(case_tumor))

## =========================
## 3) Homologous OVARY controls (GTEx_homologous)
## =========================
controls_gtex <- phenoDF %>%
  filter(
    grepl(tumor, biopsy.site, ignore.case = TRUE),
    (data.source == "GTEX" & sample.type == "normal") |
      (data.source == "TCGA" & sample.type == "adjacent")
  ) %>%
  pull(sample.id) %>%
  unique()

pheno_controls<- phenoDF %>% filter(sample.id %in% controls_gtex)
print(table(pheno_controls$data.source))

## =========================
## 4) Autoencoder-selected controls (GTEx different tissues)
## =========================
set.seed(123)
controls_autoencoder_gtex <- computeRefTissue(
  case_id      = case_tumor,
  adjacent     = FALSE,   #if TRUE = normal GTEX + adjacent TCGA & TARGET
  source       = "octad",
  control_size = 90) %>%
  as.character() %>%
  unique()

controls_autoencoder <- setdiff(controls_autoencoder_gtex, controls_gtex)
length(controls_autoencoder)

## =========================
## 5) Load counts all samples
## =========================
all_samples <- union(case_tumor, union(controls_gtex, controls_autoencoder))

counts_all <- octad::loadOctadCounts(
  sample_vector = all_samples,
  type          = "counts",
  file          = H5_PATH
)

## =========================
## 6) Get metadata
## =========================
metadata_all_raw <- phenoDF %>%
  filter(sample.id %in% colnames(counts_all)) %>%
  mutate(sample.id = as.character(sample.id))

metadata_all_raw <- metadata_all_raw[match(colnames(counts_all), metadata_all_raw$sample.id), , drop = FALSE]
stopifnot(identical(colnames(counts_all), metadata_all_raw$sample.id))

metadata_all <- metadata_all_raw %>%
  dplyr::mutate(
    group_pca = case_when(
      sample.id %in% case_tumor                ~ "TCGA_tumor",
      sample.id %in% controls_autoencoder ~ "control_autoencoder",
      sample.id %in% controls_gtex            ~ "control_GTEx"))
print(table(metadata_all$group_pca))
print(table(metadata_all$biopsy.site, metadata_all$group_pca))

## =========================
## 7) De-transform: log2(x + 1) --> raw counts
## =========================
expr_raw_all <- 2^counts_all - 1
expr_raw_all <- round(expr_raw_all)
storage.mode(expr_raw_all) <- "integer"
expr_raw_all <- as.data.frame(expr_raw_all)
## =========================
## 10) Save outputs (RDS only)
## =========================
counts_batch <- union(case_tumor, controls_gtex)
counts_auto <- union(case_tumor, controls_autoencoder)

saveRDS(expr_raw_all[counts_batch],
        file = file.path(OUT_DIR, "expr_raw_Tumor_Ctlhomolog.rds"))

saveRDS(expr_raw_all[counts_auto],
        file = file.path(OUT_DIR, "expr_raw_Tumor_CtlAuto.rds"))

saveRDS(metadata_all,
        file = file.path(OUT_DIR, "metadata_all.rds"))


