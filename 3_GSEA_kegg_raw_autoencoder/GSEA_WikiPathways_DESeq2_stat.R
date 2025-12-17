library(dplyr)
library(ggvenn)
library(ggplot2)

base_dir <- "~/Ovary_signatures/"

paths <- list(
  raw_full_deseq2  = file.path(base_dir, "0_DGE_raw",         "DE_full_OVARY_DESeq2.rds"),
  ae_full_deseq2   = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds")
)

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

DE_full_raw_DESeq2 <- read_rds_safe(paths$raw)
DE_full_ae_DESeq2  <- read_rds_safe(paths$ae)
