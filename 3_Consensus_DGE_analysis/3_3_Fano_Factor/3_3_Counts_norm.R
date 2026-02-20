suppressPackageStartupMessages({
  library(octad); library(octad.db)
  library(dplyr); library(DESeq2)
})

options(stringsAsFactors = FALSE)
tumor   <- "OVARY"
H5_PATH <- "/STORAGE/csbig/jruiz/Octad/octad.counts.and.tpm.h5"
OUT_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------- helpers --------------------
make_dds <- function(count_mat, meta_df) {
  meta_df <- meta_df[match(colnames(count_mat), meta_df$sample.id), , drop = FALSE]
  stopifnot(identical(colnames(count_mat), meta_df$sample.id))
  
  n_min <- min(table(meta_df$group_pca))
  keep  <- rowSums(count_mat >= 10) >= n_min
  count_mat <- count_mat[keep, , drop = FALSE]
  
  colData <- data.frame(condition = factor(meta_df$group_pca))
  rownames(colData) <- meta_df$sample.id
  
  DESeqDataSetFromMatrix(countData = count_mat, colData = colData, design = ~ condition)
}

# -------------------- metadata + ids --------------------
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")

case_tumor <- phenoDF %>%
  filter(grepl(tumor, biopsy.site, ignore.case = TRUE),
         data.source == "TCGA",
         sample.type != "adjacent") %>%
  pull(sample.id) %>% as.character()

controls_gtex <- phenoDF %>%
  filter(grepl(tumor, biopsy.site, ignore.case = TRUE),
         (data.source == "GTEX" & sample.type == "normal") |
           (data.source == "TCGA" & sample.type == "adjacent")) %>%
  pull(sample.id) %>% as.character() %>% unique()

set.seed(123)
controls_autoencoder_gtex <- computeRefTissue(
  case_id      = case_tumor,
  adjacent     = FALSE,
  source       = "octad",
  control_size = 90
) %>% as.character() %>% unique()

controls_autoencoder <- setdiff(controls_autoencoder_gtex, controls_gtex)

all_samples <- union(case_tumor, union(controls_gtex, controls_autoencoder))

# -------------------- load counts and metadata aligned --------------------
counts_all <- octad::loadOctadCounts(
  sample_vector = all_samples,
  type          = "counts",
  file          = H5_PATH
)

metadata_all <- phenoDF %>%
  filter(sample.id %in% colnames(counts_all)) %>%
  mutate(sample.id = as.character(sample.id),
         group_pca = case_when(
           sample.id %in% case_tumor          ~ "TCGA_tumor",
           sample.id %in% controls_autoencoder ~ "control_autoencoder",
           sample.id %in% controls_gtex       ~ "control_GTEx",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(group_pca))

metadata_all <- metadata_all[match(colnames(counts_all), metadata_all$sample.id), , drop = FALSE]
stopifnot(identical(colnames(counts_all), metadata_all$sample.id))

# -------------------- de-transform to raw counts --------------------
expr_raw_all <- round(2^counts_all - 1)
storage.mode(expr_raw_all) <- "integer"

# define sample sets (IDs)
ids_GTEX <- union(case_tumor, controls_gtex)
ids_auto <- union(case_tumor, controls_autoencoder)

# save raw + metadata (same filenames)
saveRDS(expr_raw_all[, ids_GTEX, drop = FALSE],
        file = file.path(OUT_DIR, "3_3_1_expr_raw_Tumor_Ctlhomolog.rds"))
saveRDS(expr_raw_all[, ids_auto, drop = FALSE],
        file = file.path(OUT_DIR, "3_3_1_expr_raw_Tumor_CtlAuto.rds"))
saveRDS(metadata_all,
        file = file.path(OUT_DIR, "3_3_1_metadata_all.rds"))

# -------------------- DESeq objects (shared gene universe) --------------------
meta_homolog <- metadata_all %>% filter(sample.id %in% ids_GTEX)
meta_auto    <- metadata_all %>% filter(sample.id %in% ids_auto)

dds_tumor_Ctlhomolog <- make_dds(expr_raw_all[, ids_GTEX, drop = FALSE], meta_homolog)
dds_tumor_Ctlhomolog <- DESeq(dds_tumor_Ctlhomolog)

common_genes <- rownames(dds_tumor_Ctlhomolog)  # genes kept after filtering in homolog
dds_tumor_CtlAuto <- make_dds(expr_raw_all[common_genes, ids_auto, drop = FALSE], meta_auto)
dds_tumor_CtlAuto <- DESeq(dds_tumor_CtlAuto)

# -------------------- normalized counts + save (same intent/outputs) --------------------
norm_all_ovary <- counts(dds_tumor_Ctlhomolog, normalized = TRUE)
norm_all_auto  <- counts(dds_tumor_CtlAuto,    normalized = TRUE)

saveRDS(norm_all_ovary, file.path(OUT_DIR, "3_3_2_norm_counts_tumor_GTEX.rds"))
saveRDS(norm_all_auto,  file.path(OUT_DIR, "3_3_2_norm_count_tumor_AE.rds"))

###Anotation

# --- octad.db annotation (EH7272) aligned to BOTH norm matrices ---
ga <- as.data.frame(octad.db::get_ExperimentHub_data("EH7272"))
ga$ensembl <- sub("\\..*", "", as.character(ga$ensembl))
ga <- ga[!is.na(ga$ensembl) & !duplicated(ga$ensembl), ]

gU <- union(sub("\\..*", "", rownames(norm_all_ovary)),
            sub("\\..*", "", rownames(norm_all_auto)))

annot_universe <- ga[match(gU, ga$ensembl),
                     c("ensembl","gene","type","GeneID"),
                     drop = FALSE]
rownames(annot_universe) <- gU

saveRDS(annot_universe, file.path(OUT_DIR, "3_3_2_gene_annotation_UNION_universe.rds"))
