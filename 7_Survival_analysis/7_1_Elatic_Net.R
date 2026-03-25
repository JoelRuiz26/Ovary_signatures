# ============================================================
# Elastic Net Cox – TCGA OV (OPTIMIZED FINAL)
# ============================================================

load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
  library(survival)
  library(tibble)
})

set.seed(123)

# ============================================================
# 1) LOAD GENES
# ============================================================

genes_prop_final <- readRDS("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_3_genes_prop_final.rds")
#length(genes_prop_final) #[1] 13
genes_n_final    <- readRDS("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_3_genes_n_final.rds")
#length(genes_n_final) [1] 39

# ============================================================
# 2) PREPARE MATRIX
# ============================================================

prepare_matrix <- function(genes, Count_subset, clinical_ov){
  
  genes <- intersect(genes, rownames(Count_subset))
  if(length(genes) < 2) stop("Not enough genes for model")
  
  expr_mat <- Count_subset[genes, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample.id") %>%
    mutate(patient_barcode = substr(sample.id,1,12)) %>%
    left_join(clinical_ov, by=c("patient_barcode"="bcr_patient_barcode")) %>%
    filter(!is.na(overall_survival), !is.na(deceased)) %>%
    group_by(patient_barcode) %>%
    summarise(
      across(all_of(genes), \(x) median(x, na.rm=TRUE)),
      overall_survival = overall_survival[1],
      deceased = deceased[1],
      .groups="drop"
    )
  
  x <- expr_mat %>%
    dplyr::select(all_of(genes)) %>%
    as.matrix() %>%
    scale()
  
  y <- Surv(expr_mat$overall_survival, expr_mat$deceased)
  
  list(x=x, y=y)
}

# ============================================================
# 3) ELASTIC NET FUNCTION
# ============================================================

run_elastic_net <- function(genes, Count_subset, clinical_ov, alpha_val){
  
  dat <- prepare_matrix(genes, Count_subset, clinical_ov)
  
  cvfit <- cv.glmnet(
    x = dat$x,
    y = dat$y,
    family = "cox",
    alpha = alpha_val,
    nfolds = 10
  )
  
  coef_min <- coef(cvfit, s = "lambda.min")
  genes_min <- rownames(coef_min)[coef_min[,1] != 0]
  
  list(
    model = cvfit,
    genes_min = genes_min
  )
}

# ============================================================
# 4) ALPHA OPTIMIZATION
# ============================================================

alphas <- seq(0.1, 0.9, by=0.1)

get_alpha_table <- function(genes){
  do.call(rbind, lapply(alphas, function(a){
    res <- run_elastic_net(genes, Count_subset, clinical_ov, alpha_val=a)
    
    data.frame(
      alpha = a,
      cvm = min(res$model$cvm),
      n_genes = length(res$genes_min)
    )
  }))
}

alpha_prop_df <- get_alpha_table(genes_prop_final)
alpha_n_df    <- get_alpha_table(genes_n_final)

# ============================================================
# 5) SELECT BEST ALPHA 
# ============================================================

select_best_alpha <- function(df){
  df %>%
    mutate(threshold = min(cvm) + sd(cvm)) %>%
    filter(cvm <= threshold) %>%
    arrange(n_genes) %>%
    slice(1)
}

best_prop <- select_best_alpha(alpha_prop_df)
best_n    <- select_best_alpha(alpha_n_df)

# ============================================================
# 6) FINAL MODELS
# ============================================================

enet_prop <- run_elastic_net(
  genes_prop_final,
  Count_subset,
  clinical_ov,
  alpha_val = best_prop$alpha
)

enet_n <- run_elastic_net(
  genes_n_final,
  Count_subset,
  clinical_ov,
  alpha_val = best_n$alpha
)

# ============================================================
# 7) RESULTS
# ============================================================

print(best_prop)
#alpha      cvm n_genes threshold
#1   0.2 7.163629      11  7.167514
print(enet_prop$genes_min)
#[1] "HMGB3"   "CASP2"   "FANCI"   "SPAG5"   "ZWINT"   "GMNN"    "YIF1B"  
#[8] "PAK1IP1" "CCDC167" "BLM"     "CHEK1" 

print(best_n)
#alpha      cvm n_genes threshold
#1   0.7  7.163862      25  7.166731

print(enet_n$genes_min)
#[1] "LYVE1"   "RCN2"    "HMGB3"   "MCC"     "PCDH9"   "SVEP1"   "CASP2"  
#[8] "FANCI"   "SPAG5"   "COL14A1" "CXCR4"   "ZWINT"   "SUV39H2" "HPDL"   
#[15] "NBR1"    "PDE10A"  "ANO5"    "LSAMP"   "SLC16A3" "YIF1B"   "CLDN4"  
#[22] "PRAME"   "CCDC167" "TPI1"    "CHEK1" 

save.image("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_1_Image_ElasNet.RData")


