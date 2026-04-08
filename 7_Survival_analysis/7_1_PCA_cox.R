# ============================================================
# COX MULTIVARIADO FINAL — TCGA OV (LIMPIO Y CORRECTO)
# ============================================================

load("/STORAGE/csbig/jruiz/Ovary_data/7_survival/7_0_Image_OV_GDC.RData")

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  
})

# ============================================================
# 1) EXPRESIÓN → NIVEL PACIENTE
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
# 2) CLÍNICA
# ============================================================

clinical_clean <- clinical_ov %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE)

# ============================================================
# 3) JOIN
# ============================================================

dat <- expr_patient %>%
  inner_join(clinical_clean, by = c("patient_barcode" = "bcr_patient_barcode")) %>%
  filter(!is.na(overall_survival), !is.na(deceased))

# ============================================================
# 4) VARIABLES CLÍNICAS
# ============================================================

# Residual → 3 grupos (MEJOR DEFINICIÓN SEGÚN TUS RESULTADOS)
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

# Stage (ref = Stage III)
dat <- dat %>%
  mutate(
    stage_simple = factor(as.character(stage_simple)) %>%
      droplevels() %>%
      relevel(ref = "Stage II")
  )

# Grade (ref = G1/G2)
dat <- dat %>%
  mutate(
    grade_simple = factor(as.character(grade_simple)) %>%
      droplevels() %>%
      relevel(ref = "G1/G2")
  )

# ============================================================
# 5) DATA FINAL
# ============================================================

vars_model <- c("age_dx", "stage_simple", "grade_simple", "residual_group")

dat_model <- dat %>%
  tidyr::drop_na(all_of(vars_model))

y_model <- Surv(dat_model$overall_survival, dat_model$deceased)

# ============================================================
# 6) COX MULTIVARIADO FINAL
# ============================================================

fit_final <- coxph(
  y_model ~ age_dx + stage_simple + grade_simple + residual_group,
  data = dat_model,
  x = TRUE
)

# ============================================================
# 7) RESULTADOS
# ============================================================

summary(fit_final)

# ============================================================
# 8) CHECK RÁPIDO
# ============================================================

cat("Pacientes:", nrow(dat_model), "\n")
cat("Eventos:", sum(dat_model$deceased), "\n")
cat("\nResidual groups:\n")
print(table(dat_model$residual_group))



# ============================================================
# C-INDEX CON CORRECCIÓN POR OPTIMISMO (HARRELL)
# ============================================================

library(survival)

# -------------------------------
# 1) C-index aparente
# -------------------------------
dat_model$lp <- predict(fit_final, type = "lp")

cindex_app <- concordance(
  Surv(overall_survival, deceased) ~ I(-lp),
  data = dat_model
)$concordance

# -------------------------------
# 2) Bootstrap
# -------------------------------
set.seed(123)

n_boot <- 200
optimism <- numeric(n_boot)

for(i in 1:n_boot){
  
  # sample bootstrap
  idx <- sample(1:nrow(dat_model), replace = TRUE)
  boot_data <- dat_model[idx, ]
  
  # fit en bootstrap
  fit_boot <- tryCatch(
    coxph(
      Surv(overall_survival, deceased) ~ age_dx + stage_simple + grade_simple + residual_group,
      data = boot_data
    ),
    error = function(e) NULL
  )
  
  if(is.null(fit_boot)) next
  
  # ---- performance en bootstrap (train)
  lp_boot_train <- predict(fit_boot, type = "lp")
  
  c_boot <- concordance(
    Surv(boot_data$overall_survival, boot_data$deceased) ~ I(-lp_boot_train)
  )$concordance
  
  # ---- performance en datos originales (test)
  lp_boot_test <- predict(fit_boot, newdata = dat_model, type = "lp")
  
  c_test <- concordance(
    Surv(dat_model$overall_survival, dat_model$deceased) ~ I(-lp_boot_test)
  )$concordance
  
  # optimism
  optimism[i] <- c_boot - c_test
}

optimism <- optimism[!is.na(optimism)]

# -------------------------------
# 3) Corrección
# -------------------------------
cindex_corrected <- cindex_app - mean(optimism)

cat("\nC-index aparente:\n")
print(cindex_app)

cat("\nOptimismo medio:\n")
print(mean(optimism))

cat("\nC-index corregido:\n")
print(cindex_corrected)






# ============================================================
# BLOQUE PCA
# ============================================================

# ── 1. Matriz de expresión z-score para PCA ──────────────────

genes_use <- intersect(Top_MR_DepMap, colnames(dat_model))
cat("Genes disponibles para PCA:", length(genes_use), "\n")

# Construir matriz: filas = pacientes, cols = genes z-score
mat_pca <- dat_model %>%
  dplyr::select(all_of(genes_use)) %>%
  mutate(across(everything(), ~ scale(.)[, 1]))  # z-score in-place

# Eliminar genes con varianza 0 o NA tras escalar
mat_pca <- mat_pca[, apply(mat_pca, 2, function(x) var(x, na.rm = TRUE) > 0)]
mat_pca <- mat_pca[complete.cases(mat_pca), ]

cat("Dimensiones matriz PCA:", nrow(mat_pca), "x", ncol(mat_pca), "\n")

# ── 2. PCA ───────────────────────────────────────────────────

pca_res <- prcomp(mat_pca, center = TRUE, scale. = FALSE)  
# ya está z-score → scale.=FALSE para no doble-escalar

# Varianza explicada
var_exp    <- pca_res$sdev^2
prop_var   <- var_exp / sum(var_exp)
cum_var    <- cumsum(prop_var)
n_pcs_total <- length(prop_var)

# ── 3. Elbow method (segunda derivada de varianza explicada) ──

# Método 1: Codo por segunda derivada
delta1 <- diff(prop_var)
delta2 <- diff(delta1)
elbow_idx <- which.min(delta2) + 1   # +1 por el doble diff

# Método 2: Kaiser (eigenvalue > media → equivale a var > 1/p)
kaiser_idx <- sum(var_exp > mean(var_exp))

# Método 3: Varianza acumulada ≥ 80%
var80_idx <- which(cum_var >= 0.80)[1]

cat("\n── Criterios para selección de PCs ──────────────────────\n")
cat("Elbow (2ª derivada):          PC1 a PC", elbow_idx, "\n")
cat("Kaiser (eigenvalue > media):  PC1 a PC", kaiser_idx, "\n")
cat("Varianza acumulada ≥ 80%:     PC1 a PC", var80_idx,  "\n")

# ── Decisión final (ajusta si quieres otro criterio) ──────────
n_pcs_sel <- elbow_idx   # <-- cambia a kaiser_idx o var80_idx si prefieres

cat("\n→ PCs seleccionados para Cox:", n_pcs_sel, "\n\n")


# ── 4. Plot Elbow ─────────────────────────────────────────────

df_scree <- data.frame(
  PC       = seq_len(min(20, n_pcs_total)),
  PropVar  = prop_var[1:min(20, n_pcs_total)],
  CumVar   = cum_var[1:min(20, n_pcs_total)]
)

p_scree <- ggplot(df_scree, aes(x = PC, y = PropVar * 100)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_vline(xintercept = n_pcs_sel, linetype = "dashed",
             color = "firebrick", linewidth = 1) +
  annotate("text", x = n_pcs_sel + 0.4,
           y = max(df_scree$PropVar * 100) * 0.95,
           label = paste0("Elbow = PC", n_pcs_sel),
           color = "firebrick", hjust = 0, size = 3.5) +
  theme_classic(base_size = 13) +
  labs(title = "Scree plot – PCA Top_MR_DepMap",
       x = "Principal Component", y = "% Variance Explained")

p_cum <- ggplot(df_scree, aes(x = PC, y = CumVar * 100)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 3) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = n_pcs_sel, linetype = "dashed",
             color = "firebrick", linewidth = 1) +
  theme_classic(base_size = 13) +
  labs(title = "Cumulative variance explained",
       x = "Principal Component", y = "Cumulative %")

print(p_scree + p_cum)


# ── 6. Añadir scores de PCs seleccionados a dat_model ─────────

# Scores solo para las filas de dat_model que entran en PCA
rows_pca <- as.integer(rownames(mat_pca))   # índices originales

scores_df <- as.data.frame(pca_res$x[, 1:n_pcs_sel, drop = FALSE])
colnames(scores_df) <- paste0("PC", 1:n_pcs_sel)
scores_df$row_idx   <- rows_pca

# Unir por índice de fila
dat_model$row_idx <- seq_len(nrow(dat_model))

dat_model <- dat_model %>%
  left_join(scores_df, by = "row_idx") %>%
  dplyr::select(-row_idx)

# Subset con PCs completos
dat_pca <- dat_model %>%
  tidyr::drop_na(all_of(c("age_dx","stage_simple","grade_simple",
                          "residual_group", paste0("PC", 1:n_pcs_sel))))

y_pca <- Surv(dat_pca$overall_survival, dat_pca$deceased)

cat("Pacientes en modelo PCA:", nrow(dat_pca), "\n")
cat("Eventos:", sum(dat_pca$deceased), "\n")

# ============================================================
# BLOQUE COX CON PCs
# ============================================================

# ── 7. Cox univariado por PC ──────────────────────────────────

cat("\n── Cox univariado por PC ─────────────────────────────────\n")

uni_pc <- lapply(paste0("PC", 1:n_pcs_sel), function(pc){
  
  f  <- as.formula(paste("y_pca ~", pc))
  ft <- coxph(f, data = dat_pca)
  td <- broom::tidy(ft, conf.int = TRUE, exponentiate = TRUE)
  td$pc <- pc
  td
  
}) %>% bind_rows()

print(uni_pc %>%
        dplyr::select(pc, estimate, conf.low, conf.high, p.value) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# ── 8. Seleccionar PCs significativos (p < 0.20 univariado) ───

pcs_sig <- uni_pc %>%
  filter(p.value < 0.99) %>%
  pull(pc)

if(length(pcs_sig) == 0){
  warning("Ningún PC con p<0.20 univariado. Se usan todos los seleccionados por elbow.")
  pcs_sig <- paste0("PC", 1:n_pcs_sel)
}

cat("\nPCs que entran al multivariado:", paste(pcs_sig, collapse = ", "), "\n")

# ── 9. Cox multivariado clínico + PCs ─────────────────────────

formula_pca <- as.formula(
  paste("y_pca ~ age_dx + stage_simple + grade_simple + residual_group +",
        paste(pcs_sig, collapse = " + "))
)

fit_pca <- coxph(formula_pca, data = dat_pca, x = TRUE)

cat("\n── Resumen Cox clínico + PCs ────────────────────────────\n")
print(summary(fit_pca))




# ── 5. Biplot de loadings (PC1 vs PC2) ───────────────────────

p_biplot <- fviz_pca_var(
  pca_res,
  col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE,
  title = "PCA – Gene loadings (PC1 vs PC2)"
)

print(p_biplot)


# ── Biplot PC3 vs PC4 ───────────────────────────────────────



p_biplot_34 <- fviz_pca_var(
  pca_res,
  axes = c(3, 4),
  col.var = "contrib",
  gradient.cols = c("#005F73", "#B58900", "#9B2226"),
  repel = TRUE,
  title = "PCA – Gene loadings (PC3 vs PC4)"
)

# modificar capa de labels (ggrepel)
for(i in seq_along(p_biplot_34$layers)){
  if(inherits(p_biplot_34$layers[[i]]$geom, "GeomTextRepel")){
    
    p_biplot_34$layers[[i]]$aes_params$fontface <- "bold"
    p_biplot_34$layers[[i]]$aes_params$size <- 3.5
    
    # mejorar repel (CLAVE para no sobrelapar)
    p_biplot_34$layers[[i]]$geom_params$max.overlaps <- Inf
    p_biplot_34$layers[[i]]$geom_params$box.padding <- 0.5
    p_biplot_34$layers[[i]]$geom_params$point.padding <- 0.3
  }
}

print(p_biplot_34)






# ── Loadings ────────────────────────────────────────────────

loadings <- as.data.frame(pca_res$rotation)
loadings$gene <- rownames(loadings)

# ── Contribución (%) ───────────────────────────────────────

contrib <- get_pca_var(pca_res)$contrib %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

# combinar
loadings_full <- loadings %>%
  left_join(contrib, by = "gene")


get_top_genes_pc <- function(df, pc = "PC1", n = 20){
  
  df %>%
    dplyr::select(gene, !!sym(pc)) %>%
    mutate(abs_loading = abs(.data[[pc]])) %>%
    arrange(desc(abs_loading)) %>%
    slice(1:n)
}




top_PC1 <- get_top_genes_pc(loadings_full, "PC1", 20)
top_PC2 <- get_top_genes_pc(loadings_full, "PC2", 20)
top_PC3 <- get_top_genes_pc(loadings_full, "PC3", 20)
top_PC4 <- get_top_genes_pc(loadings_full, "PC4", 20)

top_PC4


top_PC4_pos <- loadings_full %>%
  arrange(desc(PC4)) %>%
  slice(1:20)

top_PC4_neg <- loadings_full %>%
  arrange(PC4) %>%
  slice(1:20)


top_PC4_abs <- loadings_full %>%
  mutate(abs_PC4 = abs(PC4)) %>%
  arrange(desc(abs_PC4)) %>%
  slice(1:20)


library(ggplot2)

p_top_PC4 <- top_PC4_abs %>%
  mutate(gene = reorder(gene, abs_PC4)) %>%
  ggplot(aes(x = gene, y = abs_PC4, fill = PC4 > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#9B2226", "FALSE" = "#005F73"),
    labels = c("FALSE" = "Negative loading", "TRUE" = "Positive loading"),
    name = "Direction"
  ) +
  theme_classic(base_size = 13) +
  theme(
    text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    title = "Top genes contributing to PC4 (absolute loadings)",
    x = "Gene",
    y = "|Loading (PC4)|"
  )

print(p_top_PC4)





# ── 10. C-index comparativo ───────────────────────────────────

# Modelo clínico base (sobre mismo subset)
fit_clin2 <- coxph(
  y_pca ~ age_dx + stage_simple + grade_simple + residual_group,
  data = dat_pca, x = TRUE
)

lp_clin2 <- predict(fit_clin2, type = "lp")
lp_pca   <- predict(fit_pca,   type = "lp")

ci_clin <- survival::concordance(y_pca ~ I(-lp_clin2))$concordance
ci_pca  <- survival::concordance(y_pca ~ I(-lp_pca))$concordance

cat("\nC-index clínico solo:  ", round(ci_clin, 4), "\n")
cat("C-index clínico + PCs: ", round(ci_pca,  4), "\n")

# ── 11. LRT: ¿mejoran los PCs el modelo clínico? ─────────────

lrt <- anova(fit_clin2, fit_pca, test = "LRT")
cat("\nLRT clínico vs clínico + PCs:\n")
print(lrt)

# ── 12. Forest plot ───────────────────────────────────────────

forest_pca <- broom::tidy(fit_pca, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(
    term = gsub("_z$", "", term),
    
    # quitar prefijos de variables categóricas
    term = gsub("stage_simple", "", term),
    term = gsub("grade_simple", "", term),
    term = gsub("residual_group", "residual  ", term),
    
    # opcional: limpiar un poco más
    term = gsub("age_dx", "Age", term),
    
    # trim espacios
    term = trimws(term)
  )


p_forest <- ggplot(forest_pca,
                   aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high,
                      color = p.value < 0.05), size = 0.7) +
  scale_color_manual(values = c("grey60", "firebrick"),
                     labels = c("p ≥ 0.05", "p < 0.05"),
                     name = "") +
  scale_x_log10() +
  theme_classic(base_size = 13) +
  labs(title = "Multivariate Cox: Clinical + PCs (Top_MR_DepMap)",
       x = "Hazard Ratio (log scale)", y = "")

print(p_forest)






