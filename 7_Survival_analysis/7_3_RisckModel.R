# ============================================================
# Validation – Elastic Net Cox model (CLEAN)
# ============================================================

load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_1_Image_ElasNet.RData")

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(timeROC)
})

# ============================================================
# 1) EXTRACT BETAS
# ============================================================

coef_min <- coef(enet_prop$model, s = "lambda.min")

betas <- coef_min[coef_min[,1] != 0, , drop=FALSE]
genes <- rownames(betas)

betas_vec <- as.numeric(betas)
names(betas_vec) <- genes

# ============================================================
# 2) BUILD DATASET
# ============================================================

df <- data.frame(
  sample.id = colnames(Count_subset)
)

df$patient_barcode <- substr(df$sample.id, 1, 12)

# expresión
expr <- t(Count_subset[genes, df$sample.id])
expr <- scale(expr)

df <- cbind(df, expr)

# añadir clínica
df <- df %>%
  left_join(clinical_ov, by=c("patient_barcode"="bcr_patient_barcode")) %>%
  filter(!is.na(overall_survival), !is.na(deceased)) %>%
  group_by(patient_barcode) %>%
  summarise(
    across(all_of(genes), median, na.rm=TRUE),
    overall_survival = overall_survival[1],
    deceased = deceased[1],
    .groups="drop"
  )

# ============================================================
# 3) RISK SCORE
# ============================================================

X <- as.matrix(df[, genes])

df$risk_score <- as.vector(X %*% betas_vec)

# ============================================================
# 4) KAPLAN-MEIER
# ============================================================

df$risk_group <- ifelse(df$risk_score > median(df$risk_score), "High", "Low")

fit_km <- survfit(
  Surv(overall_survival, deceased) ~ risk_group,
  data = df
)

ggsurvplot(
  fit_km,
  data = df,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#D55E00", "#0072B2"),
  title = "Elastic Net Risk Score"
)

# ============================================================
# 5) COX DEL SCORE
# ============================================================

cox_model <- coxph(
  Surv(overall_survival, deceased) ~ risk_score,
  data = df
)

summary(cox_model)

# ============================================================
# 6) C-INDEX
# ============================================================

summary(cox_model)$concordance

# ============================================================
# 7) ROC (3 years)
# ============================================================

roc <- timeROC(
  T = df$overall_survival,
  delta = df$deceased,
  marker = df$risk_score,
  cause = 1,
  times = c(365*3),
  iid = TRUE
)

plot(roc, time=365*3, col="red")
roc$AUC