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

vars_model <- c("age_dx", "stage_simple", "grade_simple", "residual_group")

dat_model <- dat %>%
  tidyr::drop_na(all_of(vars_model))

y_model <- Surv(dat_model$overall_survival, dat_model$deceased)

# ============================================================
# 9) COX MULTIVARIADO CLÍNICO
# ============================================================

fit_clinical <- coxph(
  y_model ~ age_dx + stage_simple + grade_simple + residual_group,
  data = dat_model,
  x = TRUE
)

summary(fit_clinical)

####
dat_model$lp <- predict(fit_clinical, type = "lp")

library(survival)

cindex <- survival::concordance(
  y_model ~ I(-dat_model$lp)
)$concordance

cat("\nC-index:\n")
print(cindex)
#[1] 0.630805




set.seed(123)

n_boot <- 200
cindex_boot <- numeric(n_boot)

for(i in 1:n_boot){
  
  idx <- sample(1:nrow(dat_model), replace = TRUE)
  data_boot <- dat_model[idx, ]
  
  fit_boot <- tryCatch(
    coxph(
      Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group,
      data = data_boot
    ),
    error = function(e) NULL
  )
  
  if(is.null(fit_boot)) next
  
  lp_boot <- predict(fit_boot, newdata = data_boot, type = "lp")
  
  cindex_boot[i] <- survival::concordance(
    Surv(data_boot$overall_survival, data_boot$deceased) ~ I(-lp_boot)
  )$concordance
}

cindex_boot <- cindex_boot[!is.na(cindex_boot)]

cat("\nBootstrap C-index:\n")
print(mean(cindex_boot))
#[1] 0.6389215
print(sd(cindex_boot))
#[1] 0.02025426





#PCA
genes_interes <- Top_MR_DepMap

# asegurar que los genes existan
genes_use <- intersect(genes_interes, colnames(dat_model))

# escalar expresión (z-score)
for(g in genes_use){
  dat_model[[paste0(g, "_z")]] <- scale(dat_model[[g]])[,1]
}


formula_genes <- paste0(genes_use, "_z", collapse = " + ")

formula_full <- as.formula(
  paste("y_model ~ age_dx + stage_simple + grade_simple + residual_group +", formula_genes)
)

fit_combined <- coxph(
  formula_full,
  data = dat_model,
  x = TRUE
)

summary(fit_combined)
#Call:
#  coxph(formula = formula_full, data = dat_model, x = TRUE)

#n= 367, number of events= 232 

#coef exp(coef)  se(coef)      z Pr(>|z|)    
#age_dx                 0.024665  1.024971  0.006479  3.807 0.000141 ***
#  stage_simpleStage III  0.901938  2.464374  0.462392  1.951 0.051106 .  
#stage_simpleStage IV   1.107674  3.027309  0.483461  2.291 0.021956 *  
#  grade_simpleG3/G4     -0.094742  0.909608  0.204703 -0.463 0.643488    
#residual_group>20 mm   0.758652  2.135396  0.236980  3.201 0.001368 ** 
#  residual_group≤20 mm   0.593059  1.809515  0.207873  2.853 0.004331 ** 
#  HMGB3_z               -0.169417  0.844157  0.071525 -2.369 0.017854 *  
#  SUV39H2_z             -0.211092  0.809700  0.080247 -2.631 0.008525 ** 
#  BRCA2_z                0.291530  1.338474  0.074266  3.925 8.66e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#exp(coef) exp(-coef) lower .95 upper .95
#age_dx                   1.0250     0.9756    1.0120    1.0381
#stage_simpleStage III    2.4644     0.4058    0.9957    6.0995
#stage_simpleStage IV     3.0273     0.3303    1.1736    7.8087
#grade_simpleG3/G4        0.9096     1.0994    0.6090    1.3586
#residual_group>20 mm     2.1354     0.4683    1.3420    3.3978
#residual_group≤20 mm     1.8095     0.5526    1.2040    2.7196
#HMGB3_z                  0.8442     1.1846    0.7337    0.9712
#SUV39H2_z                0.8097     1.2350    0.6919    0.9476
#BRCA2_z                  1.3385     0.7471    1.1572    1.5482

#Concordance= 0.646  (se = 0.021 )
#Likelihood ratio test= 61.41  on 9 df,   p=7e-10
#Wald test            = 53.73  on 9 df,   p=2e-08
#Score (logrank) test = 55.07  on 9 df,   p=1e-08



anova(fit_clinical, fit_combined, test = "LRT")
#Analysis of Deviance Table
#Cox model: response is  y_model
#Model 1: ~ age_dx + stage_simple + grade_simple + residual_group
#Model 2: ~ age_dx + stage_simple + grade_simple + residual_group + HMGB3_z + SUV39H2_z + BRCA2_z
#loglik  Chisq Df Pr(>|Chi|)    
#1 -1127.7                         
#2 -1114.3 26.883  3   6.23e-06 ***




forest_combined <- broom::tidy(fit_combined, conf.int = TRUE, exponentiate = TRUE)

ggplot(forest_combined, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  scale_x_log10() +
  theme_classic() +
  labs(
    title = "Multivariate Cox model (clinical + genes)",
    x = "Hazard Ratio (log scale)",
    y = ""
  )




dat_model$lp_clinical  <- predict(fit_clinical, type = "lp")
dat_model$lp_combined  <- predict(fit_combined, type = "lp")

df_compare <- data.frame(
  model = c("Clinical", "Combined"),
  cindex = c(
    survival::concordance(y_model ~ I(-dat_model$lp_clinical))$concordance,
    survival::concordance(y_model ~ I(-dat_model$lp_combined))$concordance
  )
)

ggplot(df_compare, aes(x = model, y = cindex)) +
  geom_col() +
  theme_classic() +
  labs(title = "Model discrimination (C-index)", y = "C-index", x = "")


ggplot(dat_model, aes(x = lp_combined)) +
  geom_histogram(bins = 40) +
  theme_classic() +
  labs(title = "Distribution of combined risk score", x = "Linear predictor")




dat_model$risk_group <- ifelse(
  dat_model$lp_combined > median(dat_model$lp_combined),
  "High",
  "Low"
)

dat_model$risk_group <- factor(dat_model$risk_group)

library(survminer)

fit_km <- survfit(Surv(overall_survival, deceased) ~ risk_group, data = dat_model)

ggsurvplot(
  fit_km,
  data = dat_model,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  ggtheme = theme_classic()
)



library(rms)

dd <- datadist(dat_model)
options(datadist = "dd")

fit_rms <- cph(
  Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group +
    HMGB3_z + SUV39H2_z + BRCA2_z,
  data = dat_model,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

cal <- calibrate(fit_rms, method = "boot", B = 200, u = 365)

plot(cal, xlab = "Predicted survival", ylab = "Observed survival")






contrib <- broom::tidy(fit_combined) %>%
  mutate(contribution = statistic^2) %>%
  arrange(desc(contribution))

ggplot(contrib, aes(x = reorder(term, contribution), y = contribution)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = "Variable importance (Wald χ²)", x = "", y = "Contribution")




ggplot(dat_model, aes(x = lp_combined, y = overall_survival)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Risk score vs survival",
    x = "Combined risk score",
    y = "Survival time"
  )






set.seed(123)

n_boot <- 200
cindex_boot_combined <- numeric(n_boot)

for(i in 1:n_boot){
  
  idx <- sample(1:nrow(dat_model), replace = TRUE)
  data_boot <- dat_model[idx, ]
  
  fit_boot <- coxph(
    formula_full,
    data = data_boot
  )
  
  lp_boot <- predict(fit_boot, newdata = data_boot, type = "lp")
  
  cindex_boot_combined[i] <- survival::concordance(
    Surv(data_boot$overall_survival, data_boot$deceased) ~ I(-lp_boot)
  )$concordance
}

hist(cindex_boot_combined, main = "Bootstrap C-index (combined)", xlab = "C-index")
abline(v = mean(cindex_boot_combined), col = "red", lwd = 2)







library(timeROC)

# tiempos en días (ajusta si quieres 3 o 5 años)
times <- c(365, 1095, 1825)

# ROC clínico
roc_clinical <- timeROC(
  T = dat_model$overall_survival,
  delta = dat_model$deceased,
  marker = dat_model$lp_clinical,
  cause = 1,
  times = times,
  iid = TRUE
)

# ROC combinado
roc_combined <- timeROC(
  T = dat_model$overall_survival,
  delta = dat_model$deceased,
  marker = dat_model$lp_combined,
  cause = 1,
  times = times,
  iid = TRUE
)






plot(roc_clinical, time = 365, col = "blue", lwd = 2)
plot(roc_combined, time = 365, col = "red", add = TRUE, lwd = 2)

legend(
  "bottomright",
  legend = c(
    paste0("Clinical AUC = ", round(roc_clinical$AUC[1], 3)),
    paste0("Combined AUC = ", round(roc_combined$AUC[1], 3))
  ),
  col = c("blue", "red"),
  lwd = 2
)





plot(roc_clinical, time = 1095, col = "blue", lwd = 2)
plot(roc_combined, time = 1095, col = "red", add = TRUE, lwd = 2)





auc_table <- data.frame(
  Time = times,
  Clinical = roc_clinical$AUC,
  Combined = roc_combined$AUC
)

print(auc_table)







genes_eval <- c("HMGB3_z", "SUV39H2_z", "BRCA2_z")

roc_genes <- lapply(genes_eval, function(g){
  
  timeROC(
    T = dat_model$overall_survival,
    delta = dat_model$deceased,
    marker = dat_model[[g]],
    cause = 1,
    times = times,
    iid = TRUE
  )
  
})
names(roc_genes) <- genes_eval



cols <- c("darkgreen", "purple", "orange")

plot(roc_genes[[1]], time = 365, col = cols[1], lwd = 2)

for(i in 2:length(roc_genes)){
  plot(roc_genes[[i]], time = 365, col = cols[i], add = TRUE, lwd = 2)
}

legend(
  "bottomright",
  legend = paste0(
    genes_eval,
    " AUC=", round(sapply(roc_genes, function(x) x$AUC[1]), 3)
  ),
  col = cols,
  lwd = 2
)






auc_genes <- data.frame(
  Gene = genes_eval,
  AUC_1yr = sapply(roc_genes, function(x) x$AUC[1]),
  AUC_3yr = sapply(roc_genes, function(x) x$AUC[2]),
  AUC_5yr = sapply(roc_genes, function(x) x$AUC[3])
)

print(auc_genes)





dat_model$lp_combined <- predict(fit_combined, type = "lp")

cutpoint <- median(dat_model$lp_combined)

dat_model$risk_group <- ifelse(
  dat_model$lp_combined > cutpoint,
  "High risk",
  "Low risk"
)

dat_model$risk_group <- factor(dat_model$risk_group)



table(dat_model$risk_group)
table(dat_model$risk_group, dat_model$deceased)

dat_model %>%
  group_by(risk_group) %>%
  summarise(
    n = n(),
    events = sum(deceased),
    median_survival = median(overall_survival),
    mean_score = mean(lp_combined),
    .groups = "drop"
  )




library(survminer)

fit_km <- survfit(
  Surv(overall_survival, deceased) ~ risk_group,
  data = dat_model
)

ggsurvplot(
  fit_km,
  data = dat_model,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  ggtheme = theme_classic(),
  palette = c("#2E9FDF", "#E7B800"),
  title = "Kaplan-Meier: Combined model risk groups"
)







library(survminer)

cut <- surv_cutpoint(
  data = dat_model,
  time = "overall_survival",
  event = "deceased",
  variables = "lp_combined"
)

cutpoint <- cut$cutpoint$cutpoint

dat_model$risk_group_opt <- ifelse(
  dat_model$lp_combined > cutpoint,
  "High risk",
  "Low risk"
)

dat_model$risk_group_opt <- factor(dat_model$risk_group_opt)




fit_km_opt <- survfit(
  Surv(overall_survival, deceased) ~ risk_group_opt,
  data = dat_model
)

ggsurvplot(
  fit_km_opt,
  data = dat_model,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  ggtheme = theme_classic(),
  palette = c("#00AFBB", "#FC4E07"),
  title = "Kaplan-Meier: Optimized risk groups"
)


