# ============================================================
# Univariate + Multivariate Cox – TCGA OV (CORRECTO)
# ============================================================

load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(broom)
})

# ============================================================
# 1) GENES
# ============================================================

genes_prop_final <- symbol2ens$gene
cat("Genes input:", length(genes_prop_final), "\n")

genes_ok <- intersect(genes_prop_final, rownames(Count_subset))
cat("Genes tras intersect con matriz:", length(genes_ok), "\n")

# ============================================================
# 2) EXPRESIÓN → NIVEL PACIENTE
# ============================================================

expr_patient <- Count_subset[genes_ok, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample.id") %>%
  mutate(patient_barcode = substr(sample.id, 1, 12)) %>%
  group_by(patient_barcode) %>%
  summarise(
    across(all_of(genes_ok), median, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 3) CLÍNICA (1 FILA POR PACIENTE)
# ============================================================

clinical_clean <- clinical_ov %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE)

# ============================================================
# 4) JOIN FINAL
# ============================================================

dat <- expr_patient %>%
  inner_join(clinical_clean, by = c("patient_barcode" = "bcr_patient_barcode")) %>%
  filter(!is.na(overall_survival), !is.na(deceased))

# ============================================================
# 5) FACTORES CORRECTOS (SIN ordered, SIN basura)
# ============================================================

dat <- dat %>%
  mutate(
    stage_simple   = factor(as.character(stage_simple)),
    grade_simple   = factor(as.character(grade_simple)),
    residual_group = factor(as.character(residual_group))
  )

# eliminar niveles vacíos
dat <- dat %>%
  mutate(
    stage_simple   = droplevels(stage_simple),
    grade_simple   = droplevels(grade_simple),
    residual_group = droplevels(residual_group)
  )

# referencias (IMPORTANTE)
dat$stage_simple   <- relevel(dat$stage_simple, ref = "Stage II")
dat$grade_simple   <- relevel(dat$grade_simple, ref = "G1/G2")
dat$residual_group <- relevel(dat$residual_group, ref = "No residual")

cat("Pacientes:", nrow(dat), "\n")
cat("Eventos:", sum(dat$deceased), "\n")

# ============================================================
# 6) SURVIVAL
# ============================================================

y <- Surv(dat$overall_survival, dat$deceased)

# ============================================================
# 7) COX UNIVARIADO
# ============================================================

cox_results <- lapply(genes_ok, function(gene) {
  
  expr <- scale(dat[[gene]])[, 1]
  
  fit <- tryCatch(coxph(y ~ expr), error = function(e) NULL)
  
  if (is.null(fit)) {
    return(data.frame(
      gene = gene, HR = NA, HR_lo95 = NA, HR_hi95 = NA,
      z = NA, p_value = NA, n = NA, n_events = NA
    ))
  }
  
  s <- summary(fit)
  
  data.frame(
    gene = gene,
    HR = round(s$conf.int[1, "exp(coef)"], 3),
    HR_lo95 = round(s$conf.int[1, "lower .95"], 3),
    HR_hi95 = round(s$conf.int[1, "upper .95"], 3),
    z = round(s$coefficients[1, "z"], 3),
    p_value = signif(s$coefficients[1, "Pr(>|z|)"], 3),
    n = s$n,
    n_events = s$nevent
  )
})

cox_df <- bind_rows(cox_results) %>%
  mutate(p_adj = signif(p.adjust(p_value, method = "BH"), 3)) %>%
  arrange(p_value)

cox_sig <- cox_df %>% filter(p_value <= 0.05)


print(cox_sig)
