# ============================================================
# POST-ANALYSIS USING EXISTING MODELS (NO RETRAINING)
# ============================================================

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)

load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_1_Image_Elastic.RData")

valid_results <- readRDS("7_1_0_elastic_net_model_splits_only.rds")

# ============================================================
# 1) FILTRAR MODELOS CON GENES ESTABLES
# ============================================================
models_filtered <- valid_results[sapply(valid_results, function(x) {
  setequal(x$genes, stable_genes)})]

cat("Models containing all stable genes:", length(models_filtered), "\n")

#4
# ============================================================
# 2) EXTRAER MÉTRICAS (SIN REENTRENAR)
# ============================================================

results <- lapply(models_filtered, function(res) {
  
  df_test <- data.frame(
    time = res$y_test[,1],
    event = res$y_test[,2],
    risk = res$risk_test
  )
  
  # -------------------------
  # C-INDEX
  # -------------------------
  cind <- concordance(Surv(time, event) ~ risk, data = df_test)$concordance
  
  # HR + p-value
  # -------------------------
  cox_fit <- coxph(Surv(time, event) ~ risk, data = df_test)
  
  hr <- as.numeric(exp(coef(cox_fit)))
  pval <- summary(cox_fit)$coefficients[,"Pr(>|z|)"]
  
  return(data.frame(
    cindex = cind,
    hr = hr,
    pval = pval
  ))
})

results_df <- bind_rows(results)

# ============================================================
# 3) SUMMARY
# ============================================================

cat("\nC-index summary:\n")
print(summary(results_df$cindex))

cat("\nHR summary:\n")
print(summary(results_df$hr))

cat("\nP-value summary:\n")
print(summary(results_df$pval))

# ============================================================
# 4) PLOT 1: C-INDEX DISTRIBUTION
# ============================================================
res <- valid_results[[219]]
names(res)
plot(res$cvfit)
plot(res$cvfit$glmnet.fit)

p1 <- ggplot(results_df, aes(x = cindex)) +
  geom_histogram(bins = 20, fill = "#4DBBD5", color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  theme_classic(base_size = 14) +
  labs(
    title = "C-index Distribution",
    x = "C-index",
    y = "Count")

print(p1)

# ============================================================
# 5) PLOT 2: HR DISTRIBUTION
# ============================================================

p2 <- ggplot(results_df, aes(x = hr)) +
  geom_histogram(bins = 20, fill = "#E64B35", color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    title = "Hazard Ratio Distribution",
    x = "HR",
    y = "Count"
  )

print(p2)

# ============================================================
# 6) PLOT 3: VOLCANO (HR vs p-value)
# ============================================================

results_df$log_p <- -log10(results_df$pval)

p3 <- ggplot(results_df, aes(x = log2(hr), y = log_p)) +
  geom_point(alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    title = "HR vs Significance",
    x = "log2(HR)",
    y = "-log10(p-value)"
  )

print(p3)

# ============================================================
# 7) PLOT 4: BOXPLOTS
# ============================================================

p4 <- ggplot(results_df, aes(y = cindex)) +
  geom_boxplot(fill = "#4DBBD5") +
  theme_classic(base_size = 14) +
  labs(title = "C-index Variability")

print(p4)

# ============================================================
# 8) SAVE EVERYTHING
# ============================================================

saveRDS(results_df, "7_1_4_stable_models_metrics.rds")

ggsave("Cindex_distribution.pdf", p1, width = 6, height = 5)
ggsave("HR_distribution.pdf", p2, width = 6, height = 5)
ggsave("Volcano_HR.pdf", p3, width = 6, height = 5)

cat("\nAnalysis completed and saved.\n")