# ============================================================
# Elastic Net Cox – TCGA OV  |  VERSION CORREGIDA
# ============================================================
load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(survcomp)   
  library(dplyr)
  library(tibble)
  library(caret)
  library(ggplot2)
})

# ============================================================
# 1) LOAD GENES
# ============================================================
genes_prop_final <- symbol2ens$gene
cat("Genes input:", length(genes_prop_final), "\n")  # 23

# ============================================================
# 2) PREPARE DATA
# ============================================================

prepare_surv_data <- function(genes, expr_mat, clin_data) {
  
  genes <- intersect(genes, rownames(expr_mat))
  cat("Genes tras intersect con matriz:", length(genes), "\n")
  dat <- expr_mat[genes, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample.id") %>%
    mutate(patient_barcode = substr(sample.id, 1, 12)) %>%
    left_join(clin_data, by = c("patient_barcode" = "bcr_patient_barcode")) %>%
    filter(!is.na(overall_survival), !is.na(deceased)) %>%
    group_by(patient_barcode) %>%
    summarise(
      across(all_of(genes), \(x) median(x, na.rm = TRUE)),  # sintaxis moderna
      overall_survival = first(overall_survival),
      deceased         = first(deceased),
      .groups = "drop"
    )
  
  x <- dat %>% dplyr::select(all_of(genes)) %>% as.matrix()
  y <- Surv(dat$overall_survival, dat$deceased)
  
  return(list(x = x, y = y, data = dat))
}

# ============================================================
# 3) ESCALADO SIN LEAKAGE
# ============================================================

scale_train_test <- function(x_train, x_test) {
  
  mu  <- colMeans(x_train)
  sdv <- apply(x_train, 2, sd)
  sdv[sdv == 0] <- 1
  
  x_train_sc <- scale(x_train, center = mu, scale = sdv)
  x_test_sc  <- scale(x_test,  center = mu, scale = sdv)
  
  return(list(train  = x_train_sc,
              test   = x_test_sc,
              center = mu,
              scale  = sdv))
}

# ============================================================
# 4) FILTRO UNIVARIANTE
# ============================================================

univariate_filter <- function(x_train, y_train, pval_thr = 0.05) {
  
  pvals <- apply(x_train, 2, function(g) {
    tryCatch(
      summary(coxph(y_train ~ g))$coefficients[, "Pr(>|z|)"],
      error = function(e) NA_real_
    )
  })
  
  pvals <- pvals[!is.na(pvals)]
  sel   <- names(pvals[pvals < pval_thr])
  
  return(sel)
}

# ============================================================
# 5) ELASTIC NET CON SELECCIÓN DE LAMBDA ROBUSTA
#    lambda.1se si da ≥1 gen, si no → lambda.min
# ============================================================

train_elastic_net <- function(x_train, y_train, alpha = 0.5) {
  
  cvfit <- cv.glmnet(
    x_train, y_train,
    family = "cox",
    alpha  = alpha,
    nfolds = 10
  )
  
  # ── elegir lambda ──────────────────────────────────────────
  coef_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
  n_genes_1se <- sum(coef_1se[, 1] != 0)
  
  if (n_genes_1se >= 1) {
    lambda_opt  <- cvfit$lambda.1se
    lambda_rule <- "1se"
  } else {
    lambda_opt  <- cvfit$lambda.min
    lambda_rule <- "min"
  }
  # ──────────────────────────────────────────────────────────
  
  final_model <- glmnet(
    x_train, y_train,
    family = "cox",
    alpha  = alpha,
    lambda = lambda_opt
  )
  
  coef_vec       <- as.matrix(coef(final_model))
  coef_vec       <- coef_vec[coef_vec[, 1] != 0, , drop = FALSE]
  genes_selected <- rownames(coef_vec)
  
  return(list(
    model       = final_model,
    lambda      = lambda_opt,
    lambda_rule = lambda_rule,
    genes       = genes_selected,
    coefs       = coef_vec,
    cvfit       = cvfit
  ))
}

# ============================================================
# 6) PREPARAR DATOS COMPLETOS
# ============================================================

full_data   <- prepare_surv_data(genes_prop_final, Count_subset, clinical_ov)
x_full      <- full_data$x
y_full      <- full_data$y
patient_ids <- full_data$data$patient_barcode

cat("Pacientes con datos completos:", nrow(x_full), "\n")
cat("Eventos (fallecidos):", sum(y_full[, 2]), "\n")

# ============================================================
# 7) ESTRATIFICACIÓN PARA SPLITS  (CORREGIDO)
#    Cuartiles de tiempo de supervivencia entre eventos
#    Vivos = estrato 0; fallecidos se dividen en 4 cuartiles
# ============================================================

build_strata <- function(y) {
  
  strata <- rep(0L, nrow(y))           # vivos → estrato 0
  event_mask  <- y[, 2] == 1
  event_times <- y[event_mask, 1]
  
  q_cuts <- quantile(event_times, probs = c(0.25, 0.5, 0.75))
  
  strata[event_mask] <- as.integer(
    cut(event_times,
        breaks = c(-Inf, q_cuts, Inf),
        labels = FALSE)
  )
  return(strata)
}

strata_vec <- build_strata(y_full)
cat("Distribución de estratos:", table(strata_vec), "\n")

# ============================================================
# 8) LOOP PRINCIPAL  –  400 repeticiones
# ============================================================

n_reps      <- 300
all_results <- vector("list", n_reps)

for (rep in seq_len(n_reps)) {
  
  set.seed(rep)
  
  # ── 8a. Split estratificado ──────────────────────────────
  train_idx <- createDataPartition(strata_vec, p = 0.6, list = FALSE)[, 1]
  test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
  
  x_train_raw <- x_full[train_idx, ]
  y_train     <- y_full[train_idx]
  x_test_raw  <- x_full[test_idx, ]
  y_test      <- y_full[test_idx]
  
  # ── 8b. Filtro univariante SOLO en train ─────────────────
  sel_genes <- univariate_filter(x_train_raw, y_train, pval_thr = 0.999) ###set 0.05 for active filter
  
  if (length(sel_genes) < 2) {
    cat("Rep", rep, "- genes tras filtro univariante < 2, skip\n")
    next
  }
  
  x_train <- x_train_raw[, sel_genes, drop = FALSE]
  x_test  <- x_test_raw[,  sel_genes, drop = FALSE]
  
  # ── 8c. Escalado sin leakage ─────────────────────────────
  scaled <- scale_train_test(x_train, x_test)
  
  # ── 8d. Elastic Net ──────────────────────────────────────
  fit <- train_elastic_net(scaled$train, y_train, alpha = 0.5)
  
  if (length(fit$genes) == 0) {
    cat("Rep", rep, "- elastic net sin genes, skip\n")
    next
  }
  
  # ── 8e. Risk scores ──────────────────────────────────────
  #   type = "link" → log-hazard ratio lineal; mayor score = mayor riesgo
  risk_train <- as.vector(predict(fit$model, newx = scaled$train, type = "link"))
  risk_test  <- as.vector(predict(fit$model, newx = scaled$test,  type = "link"))
  
  # ── 8f. C-index en TEST (métrica principal) ───────────────
  #   concordance.index() espera: x = riesgo, tiempo, evento
  #   Nota: mayor risk → mayor riesgo → negamos para que coincida
  #   con la convención de survcomp (mayor x → peor pronóstico ✓)
  ci_obj <- tryCatch(
    concordance.index(
      x          =  risk_test,         # mayor = mayor riesgo
      surv.time  =  y_test[, 1],
      surv.event =  y_test[, 2],
      method     = "noether"
    ),
    error = function(e) NULL
  )
  
  cindex_test <- if (!is.null(ci_obj)) ci_obj$c.index else NA_real_
  cindex_se   <- if (!is.null(ci_obj)) ci_obj$se      else NA_real_
  cindex_p    <- if (!is.null(ci_obj)) ci_obj$p.value else NA_real_
  
  # ── 8g. C-index en TRAIN (para detectar sobreajuste) ─────
  ci_train <- tryCatch(
    concordance.index(
      x          =  risk_train,
      surv.time  =  y_train[, 1],
      surv.event =  y_train[, 2],
      method     = "noether"
    ),
    error = function(e) NULL
  )
  cindex_train <- if (!is.null(ci_train)) ci_train$c.index else NA_real_
  
  # ── 8h. Guardar resultados ───────────────────────────────
  all_results[[rep]] <- list(
    
    # índices y pacientes
    train_idx      = train_idx,
    test_idx       = test_idx,
    train_patients = patient_ids[train_idx],
    test_patients  = patient_ids[test_idx],
    
    # genes seleccionados en este split
    genes_univariate = sel_genes,          # tras filtro univariante
    genes_elastic    = fit$genes,          # tras elastic net
    coefs            = fit$coefs,
    
    # modelo
    lambda      = fit$lambda,
    lambda_rule = fit$lambda_rule,         # "1se" o "min"
    cvfit       = fit$cvfit,
    
    # escalado
    center = scaled$center,
    scale  = scaled$scale,
    
    # risk scores
    risk_train = risk_train,
    risk_test  = risk_test,
    
    # survival
    y_train = y_train,
    y_test  = y_test,
    
    # rendimiento
    cindex_train = cindex_train,
    cindex_test  = cindex_test,
    cindex_se    = cindex_se,
    cindex_p     = cindex_p
  )
  
  cat(sprintf(
    "Rep %3d | genes univ: %2d | elastic: %2d | lambda: %s | C-index test: %.3f\n",
    rep, length(sel_genes), length(fit$genes), fit$lambda_rule, cindex_test
  ))
}

# ============================================================
# 9) LIMPIEZA
# ============================================================

valid_results <- all_results[!sapply(all_results, is.null)]
cat("\nRepeticiones válidas:", length(valid_results), "/", n_reps, "\n")

# ============================================================
# 10) RESUMEN DE C-INDEX
# ============================================================

cindex_df <- data.frame(
  rep          = seq_along(valid_results),
  cindex_train = sapply(valid_results, `[[`, "cindex_train"),
  cindex_test  = sapply(valid_results, `[[`, "cindex_test"),
  lambda_rule  = sapply(valid_results, `[[`, "lambda_rule")
)

cat(sprintf("\nC-index TEST  → media: %.3f  SD: %.3f  mediana: %.3f\n",
            mean(cindex_df$cindex_test, na.rm = TRUE),
            sd(cindex_df$cindex_test,   na.rm = TRUE),
            median(cindex_df$cindex_test, na.rm = TRUE)))

cat(sprintf("C-index TRAIN → media: %.3f  SD: %.3f\n",
            mean(cindex_df$cindex_train, na.rm = TRUE),
            sd(cindex_df$cindex_train,   na.rm = TRUE)))

cat("\nUso de lambda:\n")
print(table(cindex_df$lambda_rule))

# ── Plot C-index distribution ─────────────────────────────
p_ci <- ggplot(cindex_df) +
  geom_density(aes(x = cindex_train, fill = "Train"), alpha = 0.4) +
  geom_density(aes(x = cindex_test,  fill = "Test"),  alpha = 0.4) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Train" = "#E64B35", "Test" = "#4DBBD5")) +
  theme_classic(base_size = 14) +
  labs(title = "C-index distribution across splits",
       x = "C-index", y = "Density", fill = NULL)

print(p_ci)
ggsave("7_1_2_cindex_distribution.pdf", p_ci, width = 6, height = 4)
ggsave("7_1_2_cindex_distribution.png", p_ci, width = 6, height = 4, dpi = 600)

# ============================================================
# 11) ESTABILIDAD DE GENES  (elastic net)
# ============================================================

genes_all  <- unlist(lapply(valid_results, `[[`, "genes_elastic"))
gene_freq  <- sort(table(genes_all) / length(valid_results), decreasing = TRUE)

gene_freq_df <- data.frame(
  gene      = names(gene_freq),
  frequency = as.numeric(gene_freq)
)

# ── Coeficiente mediano por gen ───────────────────────────
#   Para saber dirección del efecto (protector / riesgo)
coef_list <- lapply(valid_results, function(r) {
  df <- as.data.frame(r$coefs)
  df$gene <- rownames(r$coefs)
  df
})
coef_all <- bind_rows(coef_list)
colnames(coef_all)[1] <- "coef"

coef_summary <- coef_all %>%
  group_by(gene) %>%
  summarise(
    coef_median = median(coef, na.rm = TRUE),
    coef_sd     = sd(coef,     na.rm = TRUE),
    n_reps      = n(),
    .groups = "drop"
  )

gene_freq_df <- gene_freq_df %>%
  left_join(coef_summary, by = "gene") %>%
  arrange(desc(frequency))

# ── Plot estabilidad ──────────────────────────────────────
top_n_plot <- min(50, nrow(gene_freq_df))

p_st <- ggplot(gene_freq_df[1:top_n_plot, ],
               aes(x = reorder(gene, frequency),
                   y = frequency,
                   fill = coef_median > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE"  = "#E64B35",   # riesgo
                               "FALSE" = "#4DBBD5"),  # protector
                    labels = c("TRUE" = "Riesgo", "FALSE" = "Protector"),
                    name   = "Efecto") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(title = "Top Stable Genes",
       x = "Gene", y = "Selection frequency")

print(p_st)
ggsave("7_1_1_stable_genes_barplot.pdf", p_st, width = 7, height = 5)
ggsave("7_1_1_stable_genes_barplot.png", p_st, width = 7, height = 5, dpi = 600)

# ============================================================
# 12) FIRMA FINAL
# ============================================================

threshold    <- 0.3
stable_genes <- gene_freq_df %>%
  filter(frequency >= threshold) %>%
  pull(gene)

cat("\nGenes estables (freq ≥", threshold, "):\n")
print(stable_genes)

# ============================================================
# 13) GUARDADO
# ============================================================

saveRDS(
  list(
    results      = valid_results,
    gene_freq    = gene_freq_df,
    stable_genes = stable_genes,
    cindex_df    = cindex_df),
  file = "7_1_0_elastic_net_results.rds")

save.image(file = "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_1_Image_Elastic.RData")
#load(file = "/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_1_Image_Elastic.RData")
