library(vroom)
library(dplyr)

DGE_list_shared <- vroom("/home/jruiz/Ovary_signatures/2_Consensus_list/2_0_1_genes_concordant.tsv")

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
  dplyr::select(gene, n_sources,consensus_direction)

adenosine_hits

# missing genes
missing_adenosine <- setdiff(adenosine_genes_u, toupper(DGE_list_shared$gene))
missing_adenosine

# A tibble: 3 × 3
#gene    n_sources consensus_direction
#<chr>       <dbl> <chr>              
#1 ADA             2 up                 
#2 ADORA2A         2 down               
#3 ADORA3          2 up                 
#> 
#  > # missing genes
#  > missing_adenosine <- setdiff(adenosine_genes_u, toupper(DGE_list_shared$gene))
#> missing_adenosine
#[1] "ADORA1"  "ADORA2B" "SLC28A1" "SLC29A1"
#[5] "NT5E"    "ENTPD1"  "DPP4"    "ADK"    
#[9] "ADA2" 



library(dplyr)

# Load GEO RRA results
GEO <- readRDS("/STORAGE/csbig/hachepunto/adenosina/Ovary_signatures/4_DEG_GEO/GEO_ovarian_cancer_RRA_results.rds")

# Adenosine-related genes
adenosine_genes <- c(
  "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
  "SLC28A1", "SLC29A1",
  "NT5E", "ENTPD1", "DPP4",
  "ADK", "ADA", "ADA2"
)

adenosine_genes_u <- toupper(adenosine_genes)

# Genes FOUND in GEO RRA
adenosine_hits_geo <- GEO %>%
  mutate(gene_u = toupper(Gene.symbol)) %>%
  filter(gene_u %in% adenosine_genes_u) %>%
  dplyr::select(
    Gene.symbol,
    p_raw,
    p_adj,
    rank_adj,
    direction
  ) %>%
  arrange(p_raw)

adenosine_hits_geo


#     Gene.symbol      p_raw p_adj rank_adj direction
#1          ADA 0.08913409     1   7453.5        up
#2       ENTPD1 0.15853688     1   7453.5        up
#3      SLC28A1 0.23891299     1   7453.5        up
#4      SLC28A1 0.37433862     1   5219.0      down
#5         NT5E 0.45407910     1   5219.0      down
#6      SLC29A1 0.48875783     1   7453.5        up
#7      ADORA2A 0.50106032     1   7453.5        up
#8       ADORA1 0.74166441     1   5219.0      down
#9      ADORA2B 0.76477072     1   5219.0      down
#10        DPP4 0.82299735     1   5219.0      down
#11        ADA2 0.95735156     1   7453.5        up
#12      ADORA3 0.99512195     1   5219.0      down
#13      ADORA1 1.00000000     1   7453.5        up
#14         ADK 1.00000000     1   7453.5        up


