suppressPackageStartupMessages({
  library(dplyr)
  library(vroom)
  library(DESeq2)
  library(pheatmap)
  library(dynamicTreeCut)
})

IN_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_Counts_metadata"

OUT_DIR <- "~/Ovary_signatures/3_Consensus_DGE_analysis/3_2_DGE_volcano_plots"

# =========================
# 0) Load inputs
# =========================
dds_gtex <- readRDS(file.path(IN_DIR, "3_3_2_dds_tumor_Ctlhomolog.rds"))
dds_ae   <- readRDS(file.path(IN_DIR, "3_3_2_dds_tumor_CtlAuto.rds"))

core_signature <- vroom("~/Ovary_signatures/3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv") %>%
  filter(n_sources == 3)
core_genes <- unique(core_signature$ensembl)

# =========================
# Helpers
# =========================
annotation_colors = list(
  Condition = c(
    control = "#8FD694" ,
    TCGA_tumor = "#D4A373" 
  )
)

plot_heatmap_and_clusters <- function(mat, annotation_col = NULL, main_title, out_png, out_tsv) {
  dcols <- as.dist(1 - cor(mat, method = "pearson", use = "pairwise.complete.obs"))
  hcols <- hclust(dcols, method = "average")
  
  cl <- cutreeDynamic(
    dendro = hcols,
    distM  = as.matrix(dcols),
    deepSplit = 2,
    pamRespectsDendro = FALSE
  )
  ncl <- length(setdiff(unique(cl), 0))
  
  # save cluster labels
  write.table(
    data.frame(sample = names(cl), cluster = cl),
    file = out_tsv,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  # bigger + no rownames/colnames => no "black block"
  png(out_png, width = 4200, height = 2000, res = 250)
  on.exit(dev.off(), add = TRUE)
  
  pheatmap(
    mat,
    scale = "row",
    cluster_cols = hcols,
    cluster_rows = TRUE,
    annotation_col = annotation_col,
    annotation_colors = list(   # 👈 ESTA LÍNEA
      Condition = c(
        control = "#8FD694" ,
        TCGA_tumor = "#D4A373" )),
    annotation_names_col = FALSE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    border_color = NA,
    main = paste0(main_title, " | sample clusters=", ncl),
    color = colorRampPalette(c("#2C7FB8", "grey95", "#D7301F"))(100)
    #color = colorRampPalette(c("#556B2F", "#F5F5DC", "#8B1C62"))(100)
  )
  
  message(main_title, " -> identified sample clusters: ", ncl)
  invisible(list(n_clusters = ncl, sample_cluster = cl))
}


# =========================
# 1) VST from the SAVED dds
# =========================
vsd_gtex <- varianceStabilizingTransformation(dds_gtex, blind = TRUE)
vsd_ae   <- varianceStabilizingTransformation(dds_ae,   blind = TRUE)

mat_gtex_all <- assay(vsd_gtex)  # genes x samples
mat_ae_all   <- assay(vsd_ae)

# keep only your core genes (after VST)
mat_gtex_core <- mat_gtex_all[rownames(mat_gtex_all) %in% core_genes, , drop = FALSE]
mat_ae_core   <- mat_ae_all[rownames(mat_ae_all) %in% core_genes, , drop = FALSE]

# =========================
# 2) Annotation from dds (for color bar only; clustering remains unsupervised)
# =========================
anno_gtex <- data.frame(Condition = colData(dds_gtex)$condition)
rownames(anno_gtex) <- colnames(dds_gtex)
anno_gtex <- anno_gtex[colnames(mat_gtex_core), , drop = FALSE]

anno_ae <- data.frame(Condition = colData(dds_ae)$condition)
rownames(anno_ae) <- colnames(dds_ae)
anno_ae <- anno_ae[colnames(mat_ae_core), , drop = FALSE]

# ---- MODIFICACIÓN: unificar controles como "control" ----
anno_gtex$Condition <- factor(ifelse(grepl("^control", anno_gtex$Condition),
                                     "control", as.character(anno_gtex$Condition)))
anno_ae$Condition   <- factor(ifelse(grepl("^control", anno_ae$Condition),
                                     "control", as.character(anno_ae$Condition)))
# ---------------------------------------------------------

# =========================
# 3) Heatmap per dataset (TUMOR + its controls)
# =========================
res_gtex_all <- plot_heatmap_and_clusters(
  mat = mat_gtex_core,
  annotation_col = anno_gtex,
  main_title = "Core signature",
  out_png = file.path(OUT_DIR, "3_2_2_HM_core_GTEX_all.png"),
  out_tsv = file.path(OUT_DIR, "3_2_2_HM_clusters_GTEX_all.tsv")
)

res_ae_all <- plot_heatmap_and_clusters(
  mat = mat_ae_core,
  annotation_col = anno_ae,
  main_title = "Core signature",
  out_png = file.path(OUT_DIR, "3_2_2_HM_core_AE_all.png"),
  out_tsv = file.path(OUT_DIR, "3_2_2_HM_clusters_AE_all.tsv")
)


