suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

## =========================
## 0) Paths
## =========================
DDSDIR   <- "~/Ovary_signatures/3_1_DE_results/"
OUTPLOT  <- "~/Ovary_signatures/3_2_PCA_volcano/"
dir.create(OUTPLOT, showWarnings = FALSE, recursive = TRUE)

## =========================
## 1) DDS to analyze (the 4 differential analyses)
## =========================
dds_files <- list(
  auto          = "auto_dds_used.rds",
  homolog       = "homolog_dds_used.rds",
  clean_no_neg  = "clean_no_neg_dds_used.rds",
  clean_rescued = "clean_rescued_dds_used.rds"
)

## =========================
## 2) Function: Volcano plot
## =========================
make_volcano <- function(res_df, title){
  
  res_df$Type <- ifelse(res_df$log2FoldChange >= 1 & res_df$padj < 0.05, "Upregulated",
                        ifelse(res_df$log2FoldChange <= -1 & res_df$padj < 0.05, "Downregulated",
                               "Not Significant"))
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Type)) +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "gray70"
    )) +
    geom_vline(xintercept = c(-1,1), color = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
    labs(title = title, x = "log2FC", y = "-log10(padj)") +
    theme_bw(base_size = 13)
  
  return(p)
}


## =========================
## 3) Function: PCA plot
## =========================
make_pca <- function(dds, title){
  # Usar variance-stabilizing transformation (vst) en lugar de rlog
  vsd <- vst(dds, blind = FALSE)
  p <- plotPCA(vsd, intgroup = "condition") +
    ggtitle(title) +
    theme_bw(base_size = 13)
  return(p)
}


## =========================
## 4) Loop over all DDS
## =========================
summary_table <- list()

for(nm in names(dds_files)){
  
  message("\n========== Running PCA + Volcano for: ", nm, " ==========\n")
  
  dds <- readRDS(file.path(DDSDIR, dds_files[[nm]]))
  
  ## Extract DE results table
  res_file <- file.path(DDSDIR, paste0(nm, "_FULL.tsv"))
  res_df <- read.delim(res_file)
  
  ## =========================
  ## Volcano plot
  ## =========================
  volcano <- make_volcano(res_df, paste0("Volcano: ", nm))
  ggsave(
    file.path(OUTPLOT, paste0(nm, "_volcano.png")),
    volcano, width = 6, height = 5, dpi = 300
  )
  
  ## =========================
  ## PCA plot
  ## =========================
  pca <- make_pca(dds, paste0("PCA: ", nm))
  ggsave(
    file.path(OUTPLOT, paste0(nm, "_pca.png")),
    pca, width = 6, height = 5, dpi = 300
  )
  
  ## =========================
  ## Count up/down
  ## =========================
  up   <- sum(res_df$log2FoldChange >= 1  & res_df$padj < 0.05, na.rm = TRUE)
  down <- sum(res_df$log2FoldChange <= -1 & res_df$padj < 0.05, na.rm = TRUE)
  ns   <- sum(!(res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1))
  
  summary_table[[nm]] <- data.frame(
    Comparison = nm,
    Upregulated = up,
    Downregulated = down,
    NotSignificant = ns,
    TotalGenes = nrow(res_df)
  )
  
  write.table(
    summary_table[[nm]],
    file = file.path(OUTPLOT, paste0(nm, "_summary_counts.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

## =========================
## 5) Save overall summary
## =========================
summary_all <- bind_rows(summary_table)
write.table(summary_all,
            file = file.path(OUTPLOT, "ALL_DE_summaries.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

print(summary_all)
#Comparison Upregulated Downregulated NotSignificant TotalGenes
#1          auto        4919          6895          12403      24217
#2       homolog        7839          5159          11219      24217
#3  clean_no_neg        7787          4938          10954      23679
#4 clean_rescued        7839          5159          11219      24217

