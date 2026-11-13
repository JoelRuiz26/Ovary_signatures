#!/usr/bin/env Rscript

# ============================================================
# Clinical annotation of the TCGA-OV cohort for survival analysis.
#
# Downloads the BCR XML clinical supplements of TCGA-OV and builds a
# one-row-per-patient table with overall survival and the covariates
# used later in the Cox models. The GDC files are cached under
# 8_Survival_analysis/GDCdata, so the script only downloads once.
#
# Output: 8_0_clinical_TCGA_OV.rds / .tsv
# ============================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(readr)
})

options(stringsAsFactors = FALSE)

base_dir <- "~/Ovary_signatures"
out_dir  <- file.path(base_dir, "8_Survival_analysis")
gdc_dir  <- file.path(path.expand(out_dir), "GDCdata")

dir.create(path.expand(out_dir), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Clinical supplements
# ============================================================
query <- GDCquery(
  project       = "TCGA-OV",
  data.category = "Clinical",
  data.format   = "bcr xml"
)

GDCdownload(query, directory = gdc_dir)

clinical_patient  <- GDCprepare_clinic(query, "patient",   directory = gdc_dir)
clinical_followup <- GDCprepare_clinic(query, "follow_up", directory = gdc_dir)
clinical_drug     <- GDCprepare_clinic(query, "drug",      directory = gdc_dir)

# ============================================================
# 2) Overall survival
# ============================================================
# The patient table is often out of date with respect to the last
# follow-up form, so both are merged and the most recent form wins.
# The event is taken from days_to_death rather than vital_status,
# which is more consistent across TCGA releases.
last_followup <- clinical_followup %>%
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
  )

survival_tbl <- clinical_patient %>%
  dplyr::select(bcr_patient_barcode, days_to_death, days_to_last_followup) %>%
  left_join(last_followup, by = "bcr_patient_barcode") %>%
  mutate(
    days_to_death         = coalesce(days_to_death_fu, days_to_death),
    days_to_last_followup = coalesce(days_to_last_followup_fu, days_to_last_followup),
    deceased              = !is.na(days_to_death),
    overall_survival      = ifelse(deceased,
                                   as.numeric(days_to_death),
                                   as.numeric(days_to_last_followup))
  ) %>%
  filter(!is.na(overall_survival), overall_survival > 0) %>%
  dplyr::select(bcr_patient_barcode, overall_survival, deceased)

# ============================================================
# 3) Covariates
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

    stage_event_clinical_stage = trimws(as.character(stage_event_clinical_stage)),
    stage_simple = case_when(
      startsWith(stage_event_clinical_stage, "Stage IV")  ~ "Stage IV",
      startsWith(stage_event_clinical_stage, "Stage III") ~ "Stage III",
      startsWith(stage_event_clinical_stage, "Stage II")  ~ "Stage II",
      startsWith(stage_event_clinical_stage, "Stage I")   ~ "Stage I",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),

    neoplasm_histologic_grade = trimws(as.character(neoplasm_histologic_grade)),
    grade_simple = case_when(
      neoplasm_histologic_grade %in% c("G1", "G2") ~ "G1/G2",
      neoplasm_histologic_grade %in% c("G3", "G4") ~ "G3/G4",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("G1/G2", "G3/G4")),

    tumor_residual_disease = trimws(as.character(tumor_residual_disease)),
    residual_group = case_when(
      tumor_residual_disease == "No Macroscopic disease"    ~ "No residual",
      tumor_residual_disease %in% c("1-10 mm", "11-20 mm")  ~ "<=20 mm",
      tumor_residual_disease == ">20 mm"                    ~ ">20 mm",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("No residual", "<=20 mm", ">20 mm"))
  ) %>%
  dplyr::select(bcr_patient_barcode, age_dx, stage_simple, grade_simple, residual_group)

# First-line regimen, kept as two indicator variables
drug_summary <- clinical_drug %>%
  mutate(drug_name = trimws(as.character(drug_name))) %>%
  group_by(bcr_patient_barcode) %>%
  summarise(
    platinum_any = any(grepl("Carboplatin|Cisplatin", drug_name, ignore.case = TRUE), na.rm = TRUE),
    taxane_any   = any(grepl("Paclitaxel|Docetaxel|Taxol", drug_name, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    platinum_any = factor(ifelse(platinum_any, "Yes", "No"), levels = c("No", "Yes")),
    taxane_any   = factor(ifelse(taxane_any,   "Yes", "No"), levels = c("No", "Yes"))
  )

# ============================================================
# 4) Final table
# ============================================================
clinical_ov <- survival_tbl %>%
  left_join(clin_vars,    by = "bcr_patient_barcode") %>%
  left_join(drug_summary, by = "bcr_patient_barcode") %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE)

cat("[INFO] patients with overall survival:", nrow(clinical_ov),
    "| deaths:", sum(clinical_ov$deceased), "\n")
cat("[INFO] median follow-up (days):",
    round(median(clinical_ov$overall_survival)), "\n")
print(summary(clinical_ov$stage_simple))
print(summary(clinical_ov$residual_group))

saveRDS(clinical_ov, file.path(path.expand(out_dir), "8_0_clinical_TCGA_OV.rds"))
write_tsv(clinical_ov, file.path(path.expand(out_dir), "8_0_clinical_TCGA_OV.tsv"))

message("Clinical table saved to: ", out_dir)
