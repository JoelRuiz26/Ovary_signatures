#!/usr/bin/env Rscript

# ============================================================
# OBJECTIVE: Identify master regulators enriched in highly
# dependent genes (DepMap CRISPR) and generate publication-
# ready visualization
# ============================================================

suppressPackageStartupMessages({
  library(vroom)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
})

options(stringsAsFactors = FALSE)

# ============================================================
# 1) Load data
# ============================================================
###Get universe all regulones ###
names_all_regulones <- vroom("~/Ovary_signatures/5_MRA/core_mra.tsv")
names_all_regulones <- names_all_regulones %>% filter(padj_emp <= 0.05) #1174

all_regulones <- readRDS("~/Ovary_signatures/5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds")
regulons_filtered <- all_regulones[names(all_regulones) %in% names_all_regulones$TF]
length(names(regulons_filtered)) #[1] 1174
#saveRDS(names(regulons_filtered), file = "~/Ovary_signatures/6_Depmap_ovary/6_2_regulons_all.rds")


##Get significant regulones (from NES elbow analysis, top)#####
regulons <- readRDS("~/Ovary_signatures/5_MRA/5_0_Top_elbow_regulones.rds")
regulons_names <- names(regulons) #328
#saveRDS(regulons_names, file = "~/Ovary_signatures/6_Depmap_ovary/6_2_regulons_names_topNES.rds")

#Get DF of DEGs with high dependency
df <- vroom("~/Ovary_signatures/6_Depmap_ovary/6_1_0_Dependency_DEG_core.tsv",
            show_col_types = FALSE)

# ============================================================
# 2) Select highly dependent genes (DepMap threshold)
# ============================================================

df_filtered <- df %>%
  filter(median_dep <= -1) # -0.5 result in 101 genes 
                           # -1 result in 65 genes 

genes_interest <- unique(df_filtered$gene) #65

# ============================================================
# 3) Define universe (genes present in regulons)
# ============================================================

all_genes <- unique(unlist(lapply(regulons_filtered, function(x) names(x$tfmode))))
length(all_genes) #22009

# Ensure proper intersection
genes_interest <- intersect(genes_interest, all_genes) 
#65
# ============================================================
# 4) Over-representation analysis (Fisher exact test)
# ============================================================

enrich_one <- function(reg_name, regulon, genes_interest, universe) {
  
  genes_reg <- intersect(names(regulon$tfmode), universe)
  
  k <- length(intersect(genes_reg, genes_interest))  # overlap
  M <- length(genes_reg)                            # regulon size
  n <- length(genes_interest)                       # query size
  N <- length(universe)                             # universe size
  
  mat <- matrix(c(
    k,
    M - k,
    n - k,
    N - M - n + k
  ), nrow = 2)
  
  pval <- fisher.test(mat, alternative = "greater")$p.value
  
  tibble(
    regulon = reg_name,
    size_regulon = M,
    overlap = k,
    p_value = pval
  )
}

results <- imap_dfr(
  regulons_filtered,
  ~ enrich_one(.y, .x, genes_interest, all_genes))
length(results$regulon) #1174


# ============================================================
# 5) Multiple testing correction
# ============================================================

results <- results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    log10_p = -log10(p_adj)
  ) %>%
  arrange(p_adj)
#get only top NES
results <- results %>% filter(regulon %in% regulons_names)
length(results$regulon) #[1] 328

# ============================================================
# 6) Filter significant regulons
# ============================================================

results_sig <- results %>%
  filter(p_adj < 0.01, overlap >= 10)

length(results_sig$regulon) # 46 with -1 and 0.01 most extrinsec
                            # 63  "  0.5  "  0.01 
#Identify the TF DGE
TF_DGE_DepMap <- results_sig %>% filter(regulon %in% df$gene) %>% pull(regulon)
length(TF_DGE_DepMap) #23
saveRDS(TF_DGE_DepMap, file = "~/Ovary_signatures/6_Depmap_ovary/6_2_0_names_TF_DEG_TopNES.rds")

# ============================================================
# 7) Prepare data for plotting
# ============================================================

plot_df <- results_sig %>%
  mutate(
    regulon = fct_reorder(regulon, overlap),
    is_TF_DGE = regulon %in% TF_DGE_DepMap)

# ============================================================
# 8) Generate publication-ready plot
# ============================================================
label_face <- ifelse(
  levels(plot_df$regulon) %in% TF_DGE_DepMap,
  "bold",   # TF marcados → negrita
  "plain"   # resto → normal
)


p <- ggplot(plot_df, aes(x = overlap, y = regulon, fill = log10_p)) +
  
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  
  scale_fill_gradientn(
    colors = c("#deebf7", "#9ecae1", "#3182bd", "#08519c"),
    name = expression(-log[10]("adj p"))
  ) +
  
  labs(
    x = "Gene overlap",
    y = NULL,
    title = "Master regulators enriched in highly dependent genes",
    subtitle = "Master regulators enriched in highly dependent genes"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13, margin = margin(b = 12)),
    
    axis.text.y = element_text(size = 12, margin = margin(r = 2),
                               face = label_face),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 13, face = "bold"),
    
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(10, 40, 10, 5)
  ) +
  
  coord_cartesian(clip = "off")

p

# ============================================================
# 9) Get regulones interest
# ============================================================

regulons_sig <- regulons_filtered[names(regulons_filtered) %in% results_sig$regulon]


# Filter

regulons_sig_filtered <- lapply(regulons_sig, function(reg) {
  genes_reg <- names(reg$tfmode)
  genes_keep <- intersect(genes_reg, genes_interest)
  list(
    tfmode = reg$tfmode[genes_keep],
    likelihood = reg$likelihood[match(genes_keep, names(reg$tfmode))]
  )
})

# Delete empty regulones
regulons_sig_filtered <- regulons_sig_filtered[
  sapply(regulons_sig_filtered, function(x) length(x$tfmode) > 0)
]

saveRDS(
  regulons_sig_filtered,
  file = "~/Ovary_signatures/6_Depmap_ovary/6_2_1_regulons_sig_filtered.rds")


# unique vector  TFs + genes
#all_genes <- unique(c(
#  names(regulons_sig_filtered),
#  unlist(lapply(regulons_sig_filtered, function(x) names(x$tfmode)))
#))
#length(all_genes) #[1] 99

# Guardar si quieres
#saveRDS(all_genes,
#        "~/Ovary_signatures/6_Depmap_ovary/6_2_2_Important_regulones_1.rds")

# ============================================================
# 9) Save figure (high resolution, vertical layout)
# ============================================================

out_path <- "~/Ovary_signatures/6_Depmap_ovary/6_2_1_Enriched_TF_Depndency"

ggsave(paste0(out_path, ".pdf"), p,
       width = 6, height = 10, units = "in")

ggsave(paste0(out_path, ".png"), p,
       width = 6, height = 10, units = "in", dpi = 600)



