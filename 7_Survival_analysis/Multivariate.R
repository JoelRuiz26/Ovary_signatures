# ============================================================
# COX MULTIVARIADO FINAL — TCGA OV (LIMPIO Y CORRECTO)
# ============================================================ :contentReference[oaicite:0]{index=0}

# ============================================================
# 0) LIBRERÍAS Y DATOS
# ============================================================

load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggrepel)
})

# ============================================================
# 1) PROCESAMIENTO DE EXPRESIÓN (NIVEL PACIENTE)
# ============================================================

genes_ok <- intersect(symbol2ens$gene, rownames(Count_subset))

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
# 2) PROCESAMIENTO CLÍNICO
# ============================================================

clinical_clean <- clinical_ov %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE)

# ============================================================
# 3) INTEGRACIÓN EXPRESIÓN + CLÍNICA
# ============================================================

dat <- expr_patient %>%
  inner_join(clinical_clean, by = c("patient_barcode" = "bcr_patient_barcode")) %>%
  filter(!is.na(overall_survival), !is.na(deceased))

# ============================================================
# 4) CONSTRUCCIÓN DE VARIABLES CLÍNICAS
# ============================================================

## 4.1 Residual disease
dat <- dat %>%
  mutate(
    residual_group = factor(
      case_when(
        tumor_residual_disease == "No Macroscopic disease" ~ "0-No residual",
        tumor_residual_disease %in% c("1-10 mm", "11-20 mm") ~ "1-≤20mm",
        tumor_residual_disease == ">20 mm" ~ "2->20mm"
      ),
      levels = c("0-No residual", "1-≤20mm", "2->20mm")
    )
  )

## 4.2 Stage
dat <- dat %>%
  mutate(
    stage_simple = factor(as.character(stage_simple)) %>%
      droplevels() %>%
      relevel(ref = "Stage II")
  )

## 4.3 Grade
dat <- dat %>%
  mutate(
    grade_simple = factor(as.character(grade_simple)) %>%
      droplevels() %>%
      relevel(ref = "G1/G2")
  )

# ============================================================
# 5) DATA FINAL PARA MODELADO
# ============================================================

vars_model <- c("age_dx", "stage_simple", "grade_simple", "residual_group")

dat_model <- dat %>%
  tidyr::drop_na(all_of(vars_model))

y_model <- Surv(dat_model$overall_survival, dat_model$deceased)

# ============================================================
# 6) MODELO COX MULTIVARIADO
# ============================================================

fit_final <- coxph(
  y_model ~ age_dx + stage_simple + grade_simple + residual_group,
  data = dat_model,
  x = TRUE
)

# ============================================================
# 7) RESULTADOS DEL MODELO
# ============================================================

summary(fit_final)

# ============================================================
# 8) CONTROL BÁSICO DE DATOS
# ============================================================

cat("Pacientes:", nrow(dat_model), "\n")
cat("Eventos:", sum(dat_model$deceased), "\n")

cat("\nResidual groups:\n")
print(table(dat_model$residual_group))

# ============================================================
# 9) C-INDEX (CORRECCIÓN POR OPTIMISMO – HARRELL)
# ============================================================

library(survival)

## 9.1 C-index aparente
dat_model$lp <- predict(fit_final, type = "lp")

cindex_app <- concordance(
  Surv(overall_survival, deceased) ~ I(-lp),
  data = dat_model
)$concordance

## 9.2 Bootstrap
set.seed(123)

n_boot <- 200
optimism <- numeric(n_boot)

for(i in 1:n_boot){
  
  idx <- sample(1:nrow(dat_model), replace = TRUE)
  boot_data <- dat_model[idx, ]
  
  fit_boot <- tryCatch(
    coxph(
      Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group,
      data = boot_data
    ),
    error = function(e) NULL
  )
  
  if(is.null(fit_boot)) next
  
  lp_boot_train <- predict(fit_boot, type = "lp")
  
  c_boot <- concordance(
    Surv(boot_data$overall_survival, boot_data$deceased) ~ I(-lp_boot_train)
  )$concordance
  
  lp_boot_test <- predict(fit_boot, newdata = dat_model, type = "lp")
  
  c_test <- concordance(
    Surv(dat_model$overall_survival, dat_model$deceased) ~ I(-lp_boot_test)
  )$concordance
  
  optimism[i] <- c_boot - c_test
}

optimism <- optimism[!is.na(optimism)]

## 9.3 Corrección
cindex_corrected <- cindex_app - mean(optimism)

cat("\nC-index aparente:\n")
print(cindex_app)

cat("\nOptimismo medio:\n")
print(mean(optimism))

cat("\nC-index corregido:\n")
print(cindex_corrected)

# ============================================================
# 10) PCA — PREPARACIÓN
# ============================================================

genes_use <- intersect(Top_MR_DepMap, colnames(dat_model))
cat("Genes disponibles para PCA:", length(genes_use), "\n")

mat_pca <- dat_model %>%
  dplyr::select(all_of(genes_use)) %>%
  mutate(across(everything(), ~ scale(.)[, 1]))

mat_pca <- mat_pca[, apply(mat_pca, 2, function(x) var(x, na.rm = TRUE) > 0)]
mat_pca <- mat_pca[complete.cases(mat_pca), ]

cat("Dimensiones matriz PCA:", nrow(mat_pca), "x", ncol(mat_pca), "\n")

# ============================================================
# 11) PCA — CÁLCULO Y SELECCIÓN DE PCs
# ============================================================

pca_res <- prcomp(mat_pca, center = TRUE, scale. = FALSE)

var_exp  <- pca_res$sdev^2
prop_var <- var_exp / sum(var_exp)
cum_var  <- cumsum(prop_var)

delta1 <- diff(prop_var)
delta2 <- diff(delta1)
elbow_idx  <- which.min(delta2) + 1
kaiser_idx <- sum(var_exp > mean(var_exp))
var80_idx  <- which(cum_var >= 0.80)[1]

cat("\n── Criterios para selección de PCs ──────────────────────\n")
cat("Elbow:", elbow_idx, "\n")
cat("Kaiser:", kaiser_idx, "\n")
cat("Var80:", var80_idx, "\n")

n_pcs_sel <- elbow_idx

# ============================================================
# 12) PCA — INTEGRACIÓN CON DATOS CLÍNICOS
# ============================================================

rows_pca <- as.integer(rownames(mat_pca))

scores_df <- as.data.frame(pca_res$x[, 1:n_pcs_sel, drop = FALSE])
colnames(scores_df) <- paste0("PC", 1:n_pcs_sel)
scores_df$row_idx   <- rows_pca

dat_model$row_idx <- seq_len(nrow(dat_model))

dat_model <- dat_model %>%
  left_join(scores_df, by = "row_idx") %>%
  dplyr::select(-row_idx)

dat_pca <- dat_model %>%
  tidyr::drop_na(all_of(c("age_dx","stage_simple","grade_simple",
                          "residual_group", paste0("PC", 1:n_pcs_sel))))

y_pca <- Surv(dat_pca$overall_survival, dat_pca$deceased)

# ============================================================
# 13) COX CON COMPONENTES PRINCIPALES
# ============================================================

## 13.1 Univariado
uni_pc <- lapply(paste0("PC", 1:n_pcs_sel), function(pc){
  f  <- as.formula(paste("y_pca ~", pc))
  ft <- coxph(f, data = dat_pca)
  td <- broom::tidy(ft, conf.int = TRUE, exponentiate = TRUE)
  td$pc <- pc
  td
}) %>% bind_rows()

## 13.2 Selección de PCs
pcs_sig <- uni_pc %>%
  filter(p.value < 0.99) %>%
  pull(pc)

if(length(pcs_sig) == 0){
  pcs_sig <- paste0("PC", 1:n_pcs_sel)
}

## 13.3 Modelo multivariado
formula_pca <- as.formula(
  paste("y_pca ~ age_dx + stage_simple + grade_simple + residual_group +",
        paste(pcs_sig, collapse = " + "))
)

fit_pca <- coxph(formula_pca, data = dat_pca, x = TRUE)

summary(fit_pca)

# ============================================================
# 14) COMPARACIÓN DE MODELOS
# ============================================================

fit_clin2 <- coxph(
  y_pca ~ age_dx + stage_simple + grade_simple + residual_group,
  data = dat_pca, x = TRUE
)

lp_clin2 <- predict(fit_clin2, type = "lp")
lp_pca   <- predict(fit_pca,   type = "lp")

ci_clin <- survival::concordance(y_pca ~ I(-lp_clin2))$concordance
ci_pca  <- survival::concordance(y_pca ~ I(-lp_pca))$concordance

cat("\nC-index clínico solo:", ci_clin, "\n")
cat("C-index clínico + PCs:", ci_pca, "\n")

anova(fit_clin2, fit_pca, test = "LRT")






##################
# ============================================================
# PLOTS FINALES: comparación clínica vs clínica + PC4
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(broom)
  library(ggrepel)
})

# ------------------------------------------------------------
# 1) Bootstrap para comparar C-index de ambos modelos
# ------------------------------------------------------------

set.seed(123)
n_boot_cmp <- 500

boot_cmp <- lapply(seq_len(n_boot_cmp), function(b){
  
  idx <- sample(seq_len(nrow(dat_pca)), replace = TRUE)
  boot_data <- dat_pca[idx, ]
  
  fit_b_clin <- tryCatch(
    coxph(
      Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group,
      data = boot_data,
      x = TRUE
    ),
    error = function(e) NULL
  )
  
  fit_b_pca <- tryCatch(
    coxph(
      Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group + PC1 + PC2 + PC3 + PC4,
      data = boot_data,
      x = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_b_clin) || is.null(fit_b_pca)) return(NULL)
  
  lp_clin_b <- predict(fit_b_clin, newdata = dat_pca, type = "lp")
  lp_pca_b  <- predict(fit_b_pca,  newdata = dat_pca, type = "lp")
  
  cindex_clin_b <- survival::concordance(
    Surv(dat_pca$overall_survival, dat_pca$deceased) ~ I(-lp_clin_b)
  )$concordance
  
  cindex_pca_b <- survival::concordance(
    Surv(dat_pca$overall_survival, dat_pca$deceased) ~ I(-lp_pca_b)
  )$concordance
  
  tibble(
    boot = b,
    model = c("Clinical", "Clinical + PC4"),
    cindex = c(cindex_clin_b, cindex_pca_b)
  )
}) %>% bind_rows()

boot_delta <- boot_cmp %>%
  pivot_wider(names_from = model, values_from = cindex) %>%
  mutate(delta_cindex = `Clinical + PC4` - Clinical)

delta_ci <- quantile(boot_delta$delta_cindex, c(0.025, 0.5, 0.975), na.rm = TRUE)

cindex_summary <- boot_cmp %>%
  group_by(model) %>%
  summarise(
    median = median(cindex, na.rm = TRUE),
    lo = quantile(cindex, 0.025, na.rm = TRUE),
    hi = quantile(cindex, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

p_cindex_comp <- ggplot(boot_cmp, aes(x = model, y = cindex, fill = model)) +
  geom_violin(trim = FALSE, alpha = 0.25, width = 0.9, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.55, color = "grey20") +
  geom_point(
    data = cindex_summary,
    aes(x = model, y = median),
    inherit.aes = FALSE,
    shape = 21, fill = "white", color = "black", size = 3
  ) +
  geom_errorbar(
    data = cindex_summary,
    aes(x = model, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = c("Clinical" = "#005F73", "Clinical + PC4" = "#9B2226")) +
  coord_cartesian(ylim = c(min(boot_cmp$cindex, na.rm = TRUE) - 0.005,
                           max(boot_cmp$cindex, na.rm = TRUE) + 0.005)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  ) +
  labs(
    title = "Model discrimination: clinical vs clinical + PC4",
    subtitle = sprintf("Observed C-index: clinical = %.3f | clinical + PC4 = %.3f", ci_clin, ci_pca),
    x = NULL,
    y = "C-index"
  )

print(p_cindex_comp)






# ============================================================
# TIME-DEPENDENT ROC (AUC en supervivencia)
# ============================================================

library(timeROC)

# tiempos clínicamente relevantes (ajústalos si quieres)
times_eval <- c(365, 1095, 1825)  # 1, 3, 5 años

roc_clin <- timeROC(
  T = dat_pca$overall_survival,
  delta = dat_pca$deceased,
  marker = lp_clin2,
  cause = 1,
  weighting = "marginal",
  times = times_eval,
  iid = TRUE
)

roc_pca <- timeROC(
  T = dat_pca$overall_survival,
  delta = dat_pca$deceased,
  marker = lp_pca,
  cause = 1,
  weighting = "marginal",
  times = times_eval,
  iid = TRUE
)

# ------------------------------------------------------------
# Plot comparativo AUC
# ------------------------------------------------------------

df_auc <- data.frame(
  time = rep(times_eval, 2),
  AUC = c(roc_clin$AUC, roc_pca$AUC),
  model = rep(c("Clinical", "Clinical + PC4"), each = length(times_eval))
)

library(ggplot2)

p_auc <- ggplot(df_auc, aes(x = time/365, y = AUC, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#005F73", "#9B2226")) +
  theme_classic(base_size = 13) +
  labs(
    title = "Time-dependent AUC",
    subtitle = "Clinical vs Clinical + PC4",
    x = "Time (years)",
    y = "AUC"
  )

print(p_auc)




# ============================================================
# AUC GLOBAL (ROBUSTO)
# ============================================================

ci_clin_fix <- survival::concordance(
  Surv(dat_pca$overall_survival, dat_pca$deceased) ~ I(-lp_clin2),
  na.action = na.omit
)$concordance

ci_pca_fix <- survival::concordance(
  Surv(dat_pca$overall_survival, dat_pca$deceased) ~ I(-lp_pca),
  na.action = na.omit
)$concordance

df_global_auc <- data.frame(
  model = c("Clinical", "Clinical + PC4"),
  AUC = c(ci_clin_fix, ci_pca_fix)
)

print(df_global_auc)  # DEBUG CLAVE









# ============================================================
# KM basado en riesgo continuo (NO arbitrario)
# ============================================================

# eliminar NA primero (CLAVE)
df_km <- dat_pca %>%
  dplyr::filter(!is.na(lp_pca), !is.na(overall_survival), !is.na(deceased))

# definir high risk como top 25% (riesgo real alto)
cut_val <- quantile(df_km$lp_pca, 0.75, na.rm = TRUE)

df_km$risk_group <- ifelse(df_km$lp_pca >= cut_val, "High risk", "Lower risk")

table(df_km$risk_group)

# modelo KM
fit_km <- survfit(
  Surv(overall_survival, deceased) ~ risk_group,
  data = df_km
)

library(survminer)

p_km <- ggsurvplot(
  fit_km,
  data = df_km,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#005F73", "#9B2226"),
  ggtheme = theme_classic(base_size = 13),
  title = "Risk stratification based on Cox model",
  legend.title = ""
)

print(p_km)




############TRANSLATE PC4

# ============================================================
# EXTRAER LOADINGS DE PC4
# ============================================================

loadings_pc4 <- pca_res$rotation[, "PC4"]

# convertir a dataframe
df_loadings <- data.frame(
  gene = names(loadings_pc4),
  loading = loadings_pc4
)





# top 20 genes más importantes (absoluto)
top_genes_pc4 <- df_loadings %>%
  mutate(abs_loading = abs(loading)) %>%
  arrange(desc(abs_loading)) %>%
  slice(1:20)

genes_sig <- top_genes_pc4$gene



# ============================================================
# CONSTRUIR GENE SIGNATURE SCORE
# ============================================================

# matriz de expresión (z-score como en PCA)
expr_sig <- mat_pca[, genes_sig]

# vector de pesos (loadings)
weights <- top_genes_pc4$loading
names(weights) <- top_genes_pc4$gene

# asegurar mismo orden
weights <- weights[colnames(expr_sig)]

# calcular score
signature_score <- as.matrix(expr_sig) %*% weights

# añadir al dataframe
dat_pca$signature_PC4 <- as.numeric(signature_score)






cor(dat_pca$signature_PC4, dat_pca$PC4, use = "complete.obs")



fit_sig <- coxph(
  Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group + signature_PC4,
  data = dat_pca
)

summary(fit_sig)







fit_pc4 <- coxph(
  Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group + PC4,
  data = dat_pca
)

anova(fit_pc4, fit_sig, test = "LRT")



library(pheatmap)

pheatmap(
  mat_pca[, genes_sig],
  scale = "row",
  show_rownames = FALSE,
  annotation_col = data.frame(risk = dat_pca$signature_PC4)
)
