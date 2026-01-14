library(vroom)
library(dplyr)

DGE_list_shared <- vroom("/home/jruiz/Ovary_signatures/3_Consensus_DGE_analysis/3_0_1_genes_concordant.tsv")
table(DGE_list_shared$consensus_direction)
#down none   up
#2032  556 2689 


adenosine_genes <- unique(c(
  # Receptores
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  # Transportadores
  "SLC28A1",  # CNT
  "SLC29A1",  # ENT
  # Enzimas
  "NT5E",     # CD73
  "ENTPD1",   # CD39
  "DPP4",     # CD26
  "ADK",      # Adenosine kinase
  "ADA",      # Adenosine deaminase (ADA1)
  "ADA2"      # Adenosine deaminase 2
))

# define this (you were missing it)
adenosine_genes_u <- toupper(adenosine_genes)

# hits: gene + consensus_direction
adenosine_hits <- DGE_list_shared %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(gene_u %in% adenosine_genes_u) %>%
  dplyr::select(gene, n_sources, consensus_direction)

adenosine_hits
# A tibble: 4 × 3
#gene    n_sources consensus_direction
#<chr>       <dbl> <chr>              
#1 ADA             3 up                 
#2 ADORA1          2 up                 
#3 ADORA2A         2 down               
#4 ADORA3          2 up 

#gene    n_sources consensus_direction
#<chr>       <dbl> <chr>              
#1 ADA             2 up                 
#2 ADORA2A         2 down               
#3 ADORA3          2 up


# missing genes
missing_adenosine <- setdiff(adenosine_genes_u, toupper(DGE_list_shared$gene))
missing_adenosine
#[1] "ADORA2B" "SLC28A1" "SLC29A1"
#[4] "NT5E"    "ENTPD1"  "DPP4"   
#[7] "ADK"     "ADA2"


# Load GEO RRA results
GEO <- readRDS("~/Ovary_signatures/2_DEG_GEO/GEO_ovarian_cancer_RRA_results.rds")
nrow(GEO) #[1] 35764  #[1] 20656
GEO_001 <- GEO %>% filter(p_adj <= 0.01)
nrow(GEO_001) #[1] 430    #0.05 = [1] 828   

# Genes FOUND in GEO RRA
adenosine_hits_geo <- GEO %>%
  mutate(gene_u = toupper(gene)) %>%
  filter(gene_u %in% adenosine_genes_u) %>%
  dplyr::select(gene, p_raw, p_adj, direction) %>% #ADD RANKADJ
  arrange(p_raw)

adenosine_hits_geo
#Gene.symbol       p_raw     p_adj rank_adj direction
#1          ADA 0.004661855 0.1087501    821.0        up
#2       ADORA1 0.006880527 0.1410876    934.0        up
#3      SLC28A1 0.016612967 0.2517180   1264.0        up
#4       ENTPD1 0.039023343 0.4318778   1730.0        up
#5         DPP4 0.043042465 0.6969020   1026.0      down
#6         NT5E 0.184050937 1.0000000   9105.5      down
#7       ADORA1 0.314880884 1.0000000   9105.5      down
#8      SLC29A1 0.519867496 1.0000000  11251.5        up
#9      ADORA2B 1.000000000 1.0000000   9105.5      down
#10      ADORA3 1.000000000 1.0000000   9105.5      down
#11        ADA2 1.000000000 1.0000000   9105.5      down
#12     SLC28A1 1.000000000 1.0000000   9105.5      down
#13         ADK 1.000000000 1.0000000  11251.5        up
#14     ADORA2A 1.000000000 1.0000000  11251.5        up
#15        ADA2 1.000000000 1.0000000  11251.5        up
#16     ADORA2B 1.000000000 1.0000000  11251.5        up
#17      ADORA3 1.000000000 1.0000000  11251.5        up


#gene       p_raw     p_adj direction
#1      ADA 0.004661855 0.1087501        up
#2   ADORA1 0.006880527 0.1410876        up
#3  SLC28A1 0.016612967 0.2517180        up
#4   ENTPD1 0.039023343 0.4318778        up
#5     DPP4 0.043042465 0.6969020      down
#6     NT5E 0.184050937 1.0000000      down
#7  SLC29A1 0.519867496 1.0000000        up
#8      ADK 1.000000000 1.0000000        up
#9  ADORA2A 1.000000000 1.0000000        up
#10    ADA2 1.000000000 1.0000000        up
#11 ADORA2B 1.000000000 1.0000000        up
#12  ADORA3 1.000000000 1.0000000        up


###Conclusion: number ob significant genes are low, and importat genes are not included
#Lets explore each dataset
#######

GEO_full <- readRDS("~/Ovary_signatures/2_DEG_GEO/DEGs_alldsets_GEO.rds")
str(GEO_full)

# ============================================================
# GEO_full exploration (assumes: adenosine_genes_u already exists)
# ============================================================

# ---- 1) how many significant genes per dataset (unique genes only) ----
# NOTE: collapse probes per gene by keeping the probe with smallest adj.P.Val
sig_summary <- bind_rows(lapply(names(GEO_full), function(dset) {
  x <- GEO_full[[dset]] %>%
    mutate(
      gene_u = toupper(Gene.symbol),
      adj.P.Val = as.numeric(adj.P.Val)
    ) %>%
    filter(!is.na(gene_u), gene_u != "") %>%
    group_by(gene_u) %>%
    arrange(adj.P.Val, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  
  tibble(
    dataset = dset,
    n_FDR_0.05 = sum(x$adj.P.Val <= 0.05, na.rm = TRUE),
    n_FDR_0.01 = sum(x$adj.P.Val <= 0.01, na.rm = TRUE)
  )
})) %>% arrange(dataset)

sig_summary
#dataset  n_FDR_0.05 n_FDR_0.01
#<chr>         <int>      <int>
#1 GSE14407       5833       2618
#2 GSE18520      12649       9312
#3 GSE27651       9537       5137
#4 GSE38666      11917       6975
#5 GSE40595       8138       4960
#6 GSE54388       5177       2903

GEO_collapsed <- lapply(GEO_full, function(df){
  df %>%
    mutate(
      gene_u    = toupper(Gene.symbol),
      logFC     = as.numeric(logFC),
      adj.P.Val = as.numeric(adj.P.Val)
    ) %>%
    filter(!is.na(gene_u), gene_u != "") %>%
    group_by(gene_u) %>%
    arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(direction = ifelse(is.na(logFC), NA_character_, ifelse(logFC > 0, "up", "down")))
})

# ---- helper: given alpha -> intersection size + consistency table ----
consistency_by_alpha <- function(alpha){
  # genes significativos por dataset
  sig_sets <- lapply(GEO_collapsed, function(x){
    x$gene_u[!is.na(x$adj.P.Val) & x$adj.P.Val <= alpha]
  })
  
  # genes significativos en TODOS los datasets
  genes_sig_all <- Reduce(intersect, sig_sets)
  
  # dirección por dataset (solo para genes_sig_all) -> consenso + n_sources
  genes_sig_all_tbl <- bind_rows(lapply(names(GEO_collapsed), function(dset){
    GEO_collapsed[[dset]] %>%
      filter(gene_u %in% genes_sig_all, !is.na(direction)) %>%
      transmute(dataset = dset, gene_u, direction)
  })) %>%
    group_by(gene_u) %>%
    summarise(
      n_up   = sum(direction == "up"),
      n_down = sum(direction == "down"),
      consensus_direction = case_when(
        n_up > n_down ~ "up",
        n_down > n_up ~ "down",
        TRUE ~ "none"
      ),
      n_sources = pmax(n_up, n_down),
      .groups = "drop"
    ) %>%
    transmute(gene = gene_u, n_sources = as.integer(n_sources), consensus_direction) %>%
    arrange(gene)
  
  list(
    n_genes_shared_all = length(genes_sig_all),
    genes_sig_all_tbl  = genes_sig_all_tbl,
    consistency_table  = table(genes_sig_all_tbl$n_sources, genes_sig_all_tbl$consensus_direction)
  )
}

# ---- run for both cutoffs ----
res_005 <- consistency_by_alpha(0.05)
res_001 <- consistency_by_alpha(0.01)

# how many genes are shared across ALL datasets
res_005$n_genes_shared_all
#[1] 1488

res_001$n_genes_shared_all
#[1] 514

# consistency table (n_sources vs direction)
res_005$consistency_table
res_001$consistency_table
#down none  up
#3    0   32   0
#4   36    0  53
#5   72    0 130
#6  480    0 685
#> res_001$consistency_table

#down none  up
#3    0    1   0
#4    2    0   6
#5    9    0  24
#6  241    0 231


# res_005$genes_sig_all_tbl
# res_001$genes_sig_all_tbl



# ---- Adenosina dentro de los resultados (res_005 y res_001) ----
adenosine_u <- toupper(adenosine_genes)

# FDR 0.05
adenosine_in_005 <- res_005$genes_sig_all_tbl %>%
  filter(gene %in% adenosine_u)

adenosine_in_005
#gene   n_sources consensus_direction
#<chr>      <int> <chr>              
#1 ADA            6 up                 
#2 ADORA1         5 up 

# FDR 0.01
adenosine_in_001 <- res_001$genes_sig_all_tbl %>%
  filter(gene %in% adenosine_u)

adenosine_in_001
#empty


