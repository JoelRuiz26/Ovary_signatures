#!/usr/bin/env Rscript
# =============================================================================
# 5_00_expression_regulators.R
# Inputs of the regulatory network and of the Master Regulator Analysis
#
#   1. TCGA-OV expression matrix (VST) used by ARACNe-AP and aracne2regulon
#   2. List of candidate regulators given to ARACNe-AP (--tfs)
#   3. List of canonical transcription factors used to annotate regulators (isTF)
#
# Input files (data root BASE, env var OVARY_SIGNATURES):
#   data/xena/TcgaTargetGtex_gene_expected_count.gz   UCSC Xena, Toil recompute
#   data/xena/TcgaTargetGTEX_phenotype.txt.gz
#     https://toil.xenahubs.net/download/<file>
#   data/regulators/Homo_sapiens_TF.txt               AnimalTFDB, human TFs
#   data/regulators/Homo_sapiens_Cof.txt              AnimalTFDB, human cofactors
#   data/regulators/DatabaseExtract_v_1.01.txt        Lambert et al. 2018,
#     https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt
#
# Output files (BASE/5_MRA):
#   ov_tcga_vst.tsv                      genes x 427 TCGA-OV tumors
#   modulators_all_tf_tfc_mirna.txt      regulators for ARACNe-AP
#   onlyTFs.txt                          TFs with a known motif
#
# Gene annotation comes from Ensembl through biomaRt. The submitted analysis
# used the release that was current in November 2025; set ENSEMBL_VERSION to
# pin it. A different release changes a few gene symbols and biotypes.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(vroom)
  library(DESeq2)
  library(biomaRt)
})

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")
XENA    <- file.path(BASE, "data/xena")
REG_DIR <- file.path(BASE, "data/regulators")
OUT_DIR <- Sys.getenv("MRA_OUTDIR", file.path(BASE, "5_MRA"))
ENSEMBL_VERSION <- Sys.getenv("ENSEMBL_VERSION", "")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 40)


# -----------------------------------------------------------------------------
# 1. Samples: TCGA ovarian tumors and GTEx ovary
# -----------------------------------------------------------------------------

pheno <- vroom(file.path(XENA, "TcgaTargetGTEX_phenotype.txt.gz"),
               show_col_types = FALSE, .name_repair = janitor::make_clean_names)

fix_utf <- function(x) {
  x <- iconv(x, from = "", to = "UTF-8", sub = "byte")
  str_replace_all(x, "[^[:print:]]", "")
}
pheno <- pheno %>%
  mutate(across(c(study, primary_site, primary_disease_or_tissue), fix_utf)) %>%
  distinct(sample, .keep_all = TRUE)

is_tcga_ov <- !is.na(pheno$study) & pheno$study == "TCGA" &
  (str_detect(tolower(coalesce(pheno$primary_site, "")), "ovary") |
   str_detect(tolower(coalesce(pheno$primary_disease_or_tissue, "")), "ovary"))
is_gtex_ov <- !is.na(pheno$study) & pheno$study == "GTEX" &
  str_detect(tolower(coalesce(pheno$primary_site, "")), "ovary")

keep_ids <- unique(c(pheno$sample[is_tcga_ov], pheno$sample[is_gtex_ov]))


# -----------------------------------------------------------------------------
# 2. Expected counts (Toil RSEM, stored as log2(count + 1))
# -----------------------------------------------------------------------------

counts_file <- file.path(XENA, "TcgaTargetGtex_gene_expected_count.gz")
header   <- strsplit(readLines(counts_file, n = 1), "\t")[[1]]
keep_ids <- keep_ids[keep_ids %in% header]
counts   <- vroom(counts_file, col_select = all_of(c(header[1], keep_ids)),
                  show_col_types = FALSE)
exprs_matrix <- as.matrix(counts[, -1])
rownames(exprs_matrix) <- counts[[1]]
rm(counts)

pheno_df <- pheno %>%
  filter(sample %in% keep_ids) %>%
  mutate(group = factor(if_else(study == "TCGA", "Cancer", "Normal"),
                        levels = c("Normal", "Cancer"))) %>%
  arrange(match(sample, keep_ids)) %>%
  as.data.frame()
rownames(pheno_df) <- pheno_df$sample
exprs_matrix <- exprs_matrix[, pheno_df$sample, drop = FALSE]
message("Samples: ", paste(names(table(pheno_df$group)), table(pheno_df$group), collapse = " | "))
# Expected: Normal 88 | Cancer 427


# -----------------------------------------------------------------------------
# 3. Gene annotation and collapse to gene symbol
# -----------------------------------------------------------------------------

ens_clean <- sub("\\..*", "", rownames(exprs_matrix))

mart <- if (nzchar(ENSEMBL_VERSION)) {
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = as.integer(ENSEMBL_VERSION))
} else {
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
}
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype", "chromosome_name"),
  filters    = "ensembl_gene_id",
  values     = unique(ens_clean),
  mart       = mart
)

sym_vec <- setNames(gene_map$hgnc_symbol,  gene_map$ensembl_gene_id)[ens_clean]
bio_vec <- setNames(gene_map$gene_biotype, gene_map$ensembl_gene_id)[ens_clean]

keep <- !is.na(sym_vec) & sym_vec != "" & bio_vec %in% c("protein_coding", "lncRNA", "miRNA")

# NOTE: the values are summed as stored by Xena, i.e. log2(expected_count + 1),
# without back-transforming to counts. This is what the submitted analysis did.
expr_by_symbol <- rowsum(exprs_matrix[keep, , drop = FALSE], group = sym_vec[keep], reorder = FALSE)


# -----------------------------------------------------------------------------
# 4. Variance-stabilizing transformation
# -----------------------------------------------------------------------------

# VST is fitted on tumors and GTEx ovary together (blind to the design);
# only the tumor columns are written for network inference.
dds <- DESeqDataSetFromMatrix(round(expr_by_symbol), pheno_df, design = ~ group)
dds <- estimateSizeFactors(dds)
expr_vst <- assay(vst(dds, blind = TRUE))

tumor_ids <- rownames(pheno_df)[pheno_df$group == "Cancer"]
expr_vst[, tumor_ids, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  write.table(file.path(OUT_DIR, "ov_tcga_vst.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Candidate regulators for ARACNe-AP
# -----------------------------------------------------------------------------

# 5.1 miRNA genes annotated by Ensembl among the measured genes
mirna_symbols <- gene_map %>%
  filter(gene_biotype == "miRNA", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol) %>% arrange(hgnc_symbol) %>% pull(hgnc_symbol)

# 5.2 Human transcription factors (Lambert et al. 2018, "Is TF?" == Yes)
lambert <- vroom(file.path(REG_DIR, "DatabaseExtract_v_1.01.txt"),
                 .name_repair = janitor::make_clean_names, show_col_types = FALSE)
lambert_symbols <- lambert %>%
  filter(is_tf == "Yes", !is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol) %>% arrange(hgnc_symbol) %>% pull(hgnc_symbol)

# 5.3 AnimalTFDB: cofactors by symbol, TFs by Ensembl id mapped to current symbol
tfcs  <- read.delim(file.path(REG_DIR, "Homo_sapiens_Cof.txt"))$Symbol
tfs_e <- read.delim(file.path(REG_DIR, "Homo_sapiens_TF.txt"))$Ensembl
tfs_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id", values = unique(tfs_e), mart = mart)
tfs <- tfs_map$hgnc_symbol[tfs_map$hgnc_symbol != ""]

regulators_all <- unique(c(tfs, tfcs, lambert_symbols, mirna_symbols))
writeLines(regulators_all, file.path(OUT_DIR, "modulators_all_tf_tfc_mirna.txt"))
message("Candidate regulators: ", length(regulators_all))
# Expected: 4180 (4147 present in the expression matrix)


# -----------------------------------------------------------------------------
# 6. Canonical TFs with a known motif (isTF annotation in 5_0_ovary_mra.R)
# -----------------------------------------------------------------------------

only_tfs <- lambert %>%
  filter(is_tf == "Yes", tf_assessment == "Known motif") %>%
  pull(hgnc_symbol)
writeLines(only_tfs, file.path(OUT_DIR, "onlyTFs.txt"))
message("TFs with a known motif: ", length(only_tfs))
# Expected: 1107
