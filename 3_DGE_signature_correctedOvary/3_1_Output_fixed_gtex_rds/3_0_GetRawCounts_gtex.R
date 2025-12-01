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
## 0) Paths and output dir
## =========================
DIR_BASE   <- "~/Ovary_signatures/3_DGE_signature_correctedOvary/3_1_Output_fixed_gtex_rds/"
setwd(DIR_BASE)
H5_PATH    <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"

## =========================
## 1) Load OCTAD metadata
## =========================
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

## =========================
## 2) Define OVARY tumor samples (TCGA)
## =========================
case_ovary <- phenoDF %>% 
  filter(
    grepl("OVARY", biopsy.site, ignore.case = TRUE),
    data.source == "TCGA",
    sample.type != "adjacent"
  ) %>% 
  pull(sample.id)

## =========================
## 3) Homologous OVARY controls (GTEx Ovary)
## =========================
controls_ovary <- phenoDF %>%
  filter(
    grepl("OVARY", biopsy.site, ignore.case = TRUE),
    (data.source == "GTEX" & sample.type == "normal") |
      (data.source == "TCGA" & sample.type == "adjacent")
  ) %>%
  pull(sample.id) %>%
  unique()

## =========================
## 4) Autoencoder-selected controls (GTEx different tissues)
## =========================
set.seed(123)
controls_autoencoder_gtex <- computeRefTissue(
  case_id      = case_ovary,
  adjacent     = FALSE,
  source       = "octad",
  control_size = 90
) %>%
  as.character() %>%
  unique()

controls_autoencoder_gtex <- setdiff(controls_autoencoder_gtex, controls_ovary)

## =========================
## 5) Load all samples
## =========================
all_samples <- union(case_ovary, union(controls_ovary, controls_autoencoder_gtex))

counts_all <- octad::loadOctadCounts(
  sample_vector = all_samples,
  type          = "counts",
  file          = H5_PATH
)

## =========================
## 6) Align metadata
## =========================
metadata_all <- phenoDF %>%
  filter(sample.id %in% colnames(counts_all)) %>%
  mutate(sample.id = as.character(sample.id))

metadata_all <- metadata_all[match(colnames(counts_all), metadata_all$sample.id), , drop = FALSE]
stopifnot(identical(colnames(counts_all), metadata_all$sample.id))

## =========================
## 7) De-transform: log2(x + 1) --> raw counts
## =========================
expr_raw_all <- 2^counts_all - 1
storage.mode(expr_raw_all) <- "double"

## Ensure row/column order
expr_raw_all <- expr_raw_all[rownames(counts_all), colnames(counts_all)]

## =========================
## 8) Round to nearest positive integer
## =========================
## No artificial zeros, no clipping.  
expr_raw_all_int <- round(expr_raw_all)

## Use integer storage
storage.mode(expr_raw_all_int) <- "integer"

## =========================
## 9) Assign PCA groups
## =========================
metadata_all <- metadata_all %>%
  mutate(
    group_pca = case_when(
      sample.id %in% case_ovary                ~ "TCGA_tumor_ovary",
      sample.id %in% controls_autoencoder_gtex ~ "control_autoencoder_GTEx",
      sample.id %in% controls_ovary            ~ "control_ovary_GTEx",
      TRUE                                     ~ "other"
    )
  )
table(metadata_all$group_pca)
#control_autoencoder_GTEx #88 
#control_ovary_GTEx       #88 
#TCGA_tumor_ovary         #426 

table(metadata_all$biopsy.site, metadata_all$group_pca)

#                                          control_autoencoder_GTEx control_ovary_GTEx TCGA_tumor_ovary
#BRAIN - AMYGDALA                                                 4                  0                0
#BRAIN - CAUDATE (BASAL GANGLIA)                                  3                  0                0
#BRAIN - CEREBELLAR HEMISPHERE                                    4                  0                0
#BRAIN - CEREBELLUM                                              12                  0                0
#BRAIN - CORTEX                                                   1                  0                0
#BRAIN - FRONTAL CORTEX (BA9)                                     2                  0                0
#BRAIN - HIPPOCAMPUS                                              1                  0                0
#BRAIN - HYPOTHALAMUS                                             1                  0                0
#BRAIN - NUCLEUS ACCUMBENS (BASAL GANGLIA)                        4                  0                0
#BRAIN - PUTAMEN (BASAL GANGLIA)                                  1                  0                0
#BREAST - MAMMARY TISSUE                                          8                  0                0
#FALLOPIAN TUBE                                                   2                  0                0
#KIDNEY - CORTEX                                                 10                  0                0
#MUSCLE - SKELETAL                                                1                  0                0
#OVARY                                                            0                 88              426
#PITUITARY                                                       19                  0                0
#STOMACH                                                          1                  0                0
#TESTIS                                                           1                  0                0
#THYROID                                                          9                  0                0
#UTERUS                                                           4                  0                0


## =========================
## 10) Save outputs (RDS only)
## =========================
saveRDS(expr_raw_all_int,
        file = file.path(DIR_BASE, "3_0_1_expr_raw_autoGTEX.rds"))

saveRDS(metadata_all,
        file = file.path(DIR_BASE, "3_0_1_metadata_autoGTEX.rds"))


