
# ============================================================
# Survival analysis – TCGA OV
# ============================================================

setwd("~/Ovary_signatures/7_Survival_analysis/")

library(TCGAbiolinks)
library(dplyr)
library(survival)
library(survminer)
library(tibble)
library(purrr)
library(vroom)
library(ggplot2)

# ------------------------------------------------------------
# Clinical data
# ------------------------------------------------------------

query <- GDCquery(
  project="TCGA-OV",
  data.category="Clinical",
  data.format="bcr xml"
)

GDCdownload(query)

clinical_patient  <- GDCprepare_clinic(query,"patient")
clinical_followup <- GDCprepare_clinic(query,"follow_up")

clinical_ov <- clinical_patient %>%
  dplyr::select(
    bcr_patient_barcode,
    days_to_death,
    days_to_last_followup
  ) %>%
  left_join(
    clinical_followup %>%
      arrange(
        bcr_patient_barcode,
        year_of_form_completion,
        month_of_form_completion,
        day_of_form_completion
      ) %>%
      group_by(bcr_patient_barcode) %>%
      slice_tail(n=1) %>%
      ungroup() %>%
      dplyr::select(
        bcr_patient_barcode,
        days_to_death_fu = days_to_death,
        days_to_last_followup_fu = days_to_last_followup
      ),
    by="bcr_patient_barcode"
  ) %>%
  mutate(
    days_to_death = coalesce(days_to_death_fu,days_to_death),
    days_to_last_followup = coalesce(days_to_last_followup_fu,days_to_last_followup),
    deceased = !is.na(days_to_death),
    overall_survival = ifelse(
      deceased,
      as.numeric(days_to_death),
      as.numeric(days_to_last_followup)
    )
  ) %>%
  distinct(bcr_patient_barcode,.keep_all=TRUE)

# ------------------------------------------------------------
# Expression data
# ------------------------------------------------------------

Count_matrix <- readRDS(
  "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_2_norm_count_tumor_AE.rds")

metadata <- readRDS(
  "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_1_metadata_all.rds")

Top_MR_prop <- readRDS("/STORAGE/csbig/jruiz/Ovary_data/6_MRA/6_2_top5_regulons_prop_core.rds")
Top_MR_n <- readRDS("/STORAGE/csbig/jruiz/Ovary_data/6_MRA/6_1_top5_regulons_n_core.rds")

anot_genes <- readRDS("~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_2_gene_annotation_UNION_universe.rds")

metadata_tumor <- metadata %>%
  filter(!sample.type %in% c("normal","adjacent"))

Count_tumor <- Count_matrix %>%
  as.data.frame() %>%
  dplyr::select(metadata_tumor$sample.id)


# ------------------------------------------------------------
#load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")

# ------------------------------------------------------------
#Extract core genes form Regulones
# ------------------------------------------------------------
genes_prop <- unique(c(names(Top_MR_prop),unlist(lapply(Top_MR_prop, function(x) names(x$tfmode)))))
length(genes_prop) #[1] 195

genes_n    <- unique(c(names(Top_MR_n),unlist(lapply(Top_MR_n, function(x) names(x$tfmode)))))
length(genes_n) # [1] 408

genes_all <- unique(c(genes_prop, genes_n))
length(genes_all) #[1] 420  #with 32 regulones
# ============================================================
# Mapping
# ============================================================

symbol2ens <- anot_genes %>%
  dplyr::select(ensembl, gene) %>%
  filter(gene %in% genes_all)

symbol2ens <- symbol2ens %>%
  distinct(gene, .keep_all = TRUE)

# ============================================================
# Subset de matriz (rápido)
# ============================================================

Count_subset <- Count_tumor[rownames(Count_tumor) %in% symbol2ens$ensembl,]
# rename
rownames(Count_subset) <- symbol2ens$gene[match(rownames(Count_subset), symbol2ens$ensembl)]

# ============================================================
# Cox function 
# ============================================================

run_cox <- function(gene){
  
  if(!gene %in% rownames(Count_subset)) return(NULL)
  
  expr <- as.numeric(Count_subset[gene,])
  
  df <- tibble(
    sample.id = colnames(Count_subset),
    expr = expr
  ) %>%
    mutate(patient_barcode = substr(sample.id,1,12)) %>%
    left_join(
      clinical_ov,
      by=c("patient_barcode"="bcr_patient_barcode")
    ) %>%
    filter(!is.na(overall_survival), !is.na(deceased), !is.na(expr)) %>%
    group_by(patient_barcode) %>%
    summarise(
      expr = median(expr,na.rm=TRUE),
      overall_survival = overall_survival[1],
      deceased = deceased[1],
      .groups="drop"
    )
  
  if(nrow(df) < 30) return(NULL)
  if(sd(df$expr)==0) return(NULL)
  if(sum(df$deceased) < 10) return(NULL)
  
  fit <- tryCatch(
    coxph(Surv(overall_survival,deceased) ~ scale(expr), data=df),
    error=function(e) NULL
  )
  
  if(is.null(fit)) return(NULL)
  
  s <- summary(fit)
  
  tibble(
    gene = gene,
    HR = s$conf.int[,"exp(coef)"],
    CI_low = s$conf.int[,"lower .95"],
    CI_high = s$conf.int[,"upper .95"],
    pvalue = s$coefficients[,"Pr(>|z|)"]
  )
}

# ============================================================
# Cox PROP
# ============================================================

cox_prop <- map_dfr(intersect(genes_prop, rownames(Count_subset)), run_cox) %>%
  mutate(FDR = p.adjust(pvalue,"BH")) %>%
  arrange(pvalue)

vroom_write(cox_prop, "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_1_Cox_prop.tsv")

# ============================================================
# Cox N
# ============================================================

cox_n <- map_dfr(intersect(genes_n, rownames(Count_subset)), run_cox) %>%
  mutate(FDR = p.adjust(pvalue,"BH")) %>%
  arrange(pvalue)

vroom_write(cox_n, "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_2_Cox_n.tsv")


###Filter predicted genes ####
cox_n_filt <- cox_n %>%
  filter(
    pvalue < 0.05,
#    (CI_low > 1 | CI_high < 1),
#    (HR > 1.2 | HR < 0.8)
  ) %>%
  arrange(pvalue)

cox_prop_filt <- cox_prop %>%
  filter(
    pvalue < 0.05,
#    (CI_low > 1 | CI_high < 1),
#    (HR > 1.2 | HR < 0.8)
  ) %>%
  arrange(pvalue)

####

genes_prop_final <- cox_prop_filt$gene
genes_n_final    <- cox_n_filt$gene

saveRDS(genes_prop_final, "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_3_genes_prop_final.rds")
saveRDS(genes_n_final,    "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_3_genes_n_final.rds")

#save.image("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")
