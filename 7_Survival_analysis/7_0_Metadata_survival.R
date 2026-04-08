# ============================================================
# Clinical data 
# ============================================================

library(TCGAbiolinks)
library(dplyr)

# ------------------------------------------------------------
# Query
# ------------------------------------------------------------
query <- GDCquery(
  project       = "TCGA-OV",
  data.category = "Clinical",
  data.format   = "bcr xml"
)

GDCdownload(query)

clinical_patient  <- GDCprepare_clinic(query, "patient")
clinical_followup <- GDCprepare_clinic(query, "follow_up")
clinical_drug     <- GDCprepare_clinic(query, "drug")

# ============================================================
# 1 — SURVIVAL (ROBUST TCGA)
# ============================================================

clinical_ov_base <- clinical_patient %>%
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
      slice_tail(n = 1) %>%
      ungroup() %>%
      dplyr::select(
        bcr_patient_barcode,
        days_to_death_fu         = days_to_death,
        days_to_last_followup_fu = days_to_last_followup
      ),
    by = "bcr_patient_barcode"
  ) %>%
  mutate(
    days_to_death         = coalesce(days_to_death_fu, days_to_death),
    days_to_last_followup = coalesce(days_to_last_followup_fu, days_to_last_followup),
    
    # Evento (más robusto que vital_status)
    deceased = !is.na(days_to_death),
    
    overall_survival = case_when(
      deceased ~ as.numeric(days_to_death),
      !deceased ~ as.numeric(days_to_last_followup)
    )
  ) %>% filter(!is.na(overall_survival))

# ============================================================
# 2 — CLINICAL VARIABLES (SOLO LAS IMPORTANTES)
# ============================================================
clin_vars <- clinical_patient %>%
  dplyr::select(
    bcr_patient_barcode,
    age_dx = age_at_initial_pathologic_diagnosis,
    stage_event_clinical_stage,
    neoplasm_histologic_grade,
    tumor_residual_disease
  ) %>%
  mutate(
    age_dx = as.numeric(age_dx),
    
    # -----------------------------
    # STAGE
    # -----------------------------
    stage_event_clinical_stage = trimws(as.character(stage_event_clinical_stage)),
    
    stage_simple = case_when(
      startsWith(stage_event_clinical_stage, "Stage IV")  ~ "Stage IV",
      startsWith(stage_event_clinical_stage, "Stage III") ~ "Stage III",
      startsWith(stage_event_clinical_stage, "Stage II")  ~ "Stage II",
      startsWith(stage_event_clinical_stage, "Stage I")   ~ "Stage I",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("Stage I","Stage II","Stage III","Stage IV"), ordered = TRUE),
    
    # -----------------------------
    # GRADE
    # -----------------------------
    neoplasm_histologic_grade = trimws(as.character(neoplasm_histologic_grade)),
    
    grade_simple = case_when(
      neoplasm_histologic_grade %in% c("G1","G2") ~ "G1/G2",
      neoplasm_histologic_grade %in% c("G3","G4") ~ "G3/G4",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("G1/G2","G3/G4")),
    
    # -----------------------------
    # RESIDUAL
    # -----------------------------
    tumor_residual_disease = trimws(as.character(tumor_residual_disease)),
    
    residual_group = case_when(
      tumor_residual_disease == "No Macroscopic disease" ~ "No residual",
      tumor_residual_disease %in% c("1-10 mm","11-20 mm") ~ "≤20 mm",
      tumor_residual_disease == ">20 mm" ~ ">20 mm",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("No residual","≤20 mm",">20 mm"))
  )

# ============================================================
# 3 — TRATAMIENTO (MINIMO CLÍNICAMENTE INTERPRETABLE)
# ============================================================

clinical_drug <- clinical_drug %>%
  mutate(
    drug_name = trimws(as.character(drug_name)))

drug_summary <- clinical_drug %>%
  group_by(bcr_patient_barcode) %>%
  summarise(
    platinum_any = any(grepl("Carboplatin|Cisplatin", drug_name, ignore.case = TRUE), na.rm = TRUE),
    taxane_any   = any(grepl("Paclitaxel|Docetaxel|Taxol", drug_name, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    platinum_any = factor(ifelse(platinum_any, "Yes","No"), levels = c("No","Yes")),
    taxane_any   = factor(ifelse(taxane_any, "Yes","No"), levels = c("No","Yes"))
  )

# ============================================================
# 4 — FINAL DATASET PARA COX
# ============================================================

clinical_ov <- clinical_ov_base %>%
  left_join(clin_vars, by = "bcr_patient_barcode") %>%
  left_join(drug_summary, by = "bcr_patient_barcode") %>%
  filter(
    !is.na(stage_simple),
    !is.na(residual_group)
  )
# ------------------------------------------------------------
# Expression data
# ------------------------------------------------------------

Count_matrix <- readRDS(
  "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_2_norm_counts_tumor_ONLY.rds")
#[1] 60498   426
metadata <- readRDS(
  "~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_1_metadata_all.rds")


#Choose only one

#MR_high dependency MR tf
#MR_DepMap <- readRDS("~/Ovary_signatures/6_Depmap_ovary/6_4_regulons_sig_filtered.rds")
#Top_MR_DepMap <- names(MR_DepMap) 
#length(Top_MR_DepMap) #46 tf

#Genes core in high_dependency MR
#Top_MR_DepMap <- readRDS("~/Ovary_signatures/6_Depmap_ovary/6_5_Important_regulones_0.5.rds")
#length(Top_MR_DepMap) #151

#All names of TF in regulones
Top_MR_DepMap <- readRDS("~/Ovary_signatures/6_Depmap_ovary/6_4_regulons_sig_filtered.rds")
Top_MR_DepMap <- names(Top_MR_DepMap)
length(Top_MR_DepMap) #46

#Top_MR_DepMap <- readRDS(file = "~/Ovary_signatures/6_Depmap_ovary/6_3_1_names_TF_DEG_TopNES.rds")
#length(Top_MR_DepMap) #23


#All Genes core signature
#Top_MR_DepMap <- vroom("~/Ovary_signatures/3_Consensus_DGE_analysis/3_0_1_genes_concordant.tsv")
#Top_MR_DepMap <- Top_MR_DepMap %>% filter(n_sources == 3)
#length(Top_MR_DepMap$gene) #648


anot_genes <- readRDS("~/Ovary_signatures/3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_2_gene_annotation_UNION_universe.rds")

metadata_tumor <- metadata %>%
  filter(!sample.type %in% c("normal","adjacent"))

Count_tumor <- Count_matrix %>%
  as.data.frame() %>%
  dplyr::select(metadata_tumor$sample.id)

# ============================================================
# Mapping
# ============================================================

symbol2ens <- anot_genes %>%
  dplyr::select(ensembl, gene) %>%
  filter(gene %in% Top_MR_DepMap)
#  filter(gene %in% Top_MR_DepMap$gene)

symbol2ens <- symbol2ens %>%
  distinct(gene, .keep_all = TRUE)

# ============================================================
# Subset de matriz
# ============================================================

Count_subset <- Count_tumor[rownames(Count_tumor) %in% symbol2ens$ensembl,]
dim(Count_subset) #[1] 45 426

# rename
rownames(Count_subset) <- symbol2ens$gene[match(rownames(Count_subset), symbol2ens$ensembl)]

saveRDS(symbol2ens$gene, "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_3_genes_DepMap_MR.rds")

# ------------------------------------------------------------
save.image("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")


#load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")
