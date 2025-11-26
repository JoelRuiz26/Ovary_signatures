###############################################
## Matriz de eigen-espacio y Matriz_ruido
## - PCA separado en:
##     (a) OVARY tumores + controles homólogos
##     (b) OVARY tumores + controles autoencoder
## - Matriz_ruido = eigenvectores_homolog - eigenvectores_autoenc
###############################################

suppressPackageStartupMessages({
  library(stats)  # prcomp
})

options(stringsAsFactors = FALSE)

## =========================
## 0) Paths
## =========================
DIR_BASE <- "~/Ovary_signatures/3_DGE_signature_correctedOvary/"
OUT_DIR  <- file.path(DIR_BASE, "3_1_Output_fixed_gtex_rds/")

## =========================
## 1) Cargar matrices de conteos filtrados
## =========================
expr_raw_ovary_homolog <- readRDS(
  file.path(OUT_DIR, "Ovary_TCGA_tumor_plus_homologousControls_rawCounts_filtered.rds")
)

expr_raw_ovary_autoenc <- readRDS(
  file.path(OUT_DIR, "Ovary_TCGA_tumor_plus_autoencoderControls_rawCounts_filtered.rds")
)

## Asegurar mismos genes (mismas filas y mismo orden)
common_genes <- intersect(
  rownames(expr_raw_ovary_homolog),
  rownames(expr_raw_ovary_autoenc)
)

expr_raw_ovary_homolog <- expr_raw_ovary_homolog[common_genes, , drop = FALSE]
expr_raw_ovary_autoenc <- expr_raw_ovary_autoenc[common_genes, , drop = FALSE]

## =========================
## 2) PCA por separado en cada matriz
##    (samples en filas => transponer)
## =========================
pca_homolog <- prcomp(
  t(expr_raw_ovary_homolog),
  center = TRUE,
  scale. = TRUE
)

pca_autoenc <- prcomp(
  t(expr_raw_ovary_autoenc),
  center = TRUE,
  scale. = TRUE
)

## =========================
## 3) Extraer eigenvectores (loadings / rotation)
## =========================
loadings_homolog <- pca_homolog$rotation   # genes x PCs
loadings_autoenc <- pca_autoenc$rotation   # genes x PCs

## Asegurar mismo conjunto de genes (debería coincidir con common_genes)
common_genes_rot <- intersect(
  rownames(loadings_homolog),
  rownames(loadings_autoenc)
)
stopifnot(identical(common_genes_rot, common_genes))
loadings_homolog <- loadings_homolog[common_genes_rot, , drop = FALSE]
loadings_autoenc <- loadings_autoenc[common_genes_rot, , drop = FALSE]

## =========================
## 4) Alinear número de componentes (usar TODOS los PCs)
## =========================

## Verificar que ambos PCA tienen exactamente el mismo número de PCs
stopifnot(ncol(loadings_homolog) == ncol(loadings_autoenc))

## Por claridad, usamos un nombre k pero es solo n_PCs
k <- ncol(loadings_homolog)

## No recortamos nada: usamos todas las columnas tal cual
loadings_homolog_k <- loadings_homolog
loadings_autoenc_k <- loadings_autoenc

## =========================
## 5) Matriz_ruido = eigenvectores_homolog - eigenvectores_autoenc
## =========================
Matriz_ruido <- loadings_homolog_k - loadings_autoenc_k

## =========================
## 6) Guardar eigen-espacio original y Matriz_ruido
## =========================
saveRDS(
  loadings_homolog_k,
  file = file.path(OUT_DIR, "Ovary_eigenspace_homolog_loadings_firstKPCs.rds")
)

saveRDS(
  loadings_autoenc_k,
  file = file.path(OUT_DIR, "Ovary_eigenspace_autoenc_loadings_firstKPCs.rds")
)

saveRDS(
  Matriz_ruido,
  file = file.path(OUT_DIR, "Ovary_Matriz_ruido_PCA_loadings_homolog_minus_autoenc.rds")
)

## =========================
## 7) Matriz homóloga corregida en eigen-espacio
##    (loadings corregidos con tu fórmula)
## =========================

# Asegurarse de que Matriz_ruido y loadings_homolog_k tienen misma dim
stopifnot(
  all(dim(Matriz_ruido) == dim(loadings_homolog_k)),
  all(rownames(Matriz_ruido) == rownames(loadings_homolog_k))
)

# Aplicar: loadings_homolog_corregidos = sqrt( loadings_homolog_k^2 - Matriz_ruido^2 )
loadings_homolog_corregidos <- sqrt(
  pmax(loadings_homolog_k^2 - Matriz_ruido^2, 0)
)

## =========================
## 8) Diagnóstico de eigenvectores unusual
## =========================

# 1) NaN / Inf / imaginarios
num_nan <- sum(is.nan(loadings_homolog_corregidos))
num_inf <- sum(is.infinite(loadings_homolog_corregidos))
tiene_parte_imaginaria <- any(Im(loadings_homolog_corregidos) != 0)

# 2) Rango de valores (en PCA con datos centrados y escalados
#    los loadings originales suelen estar entre -1 y 1, aquí solo >= 0)
min_val  <- min(loadings_homolog_corregidos, na.rm = TRUE)
max_val  <- max(loadings_homolog_corregidos, na.rm = TRUE)

cat("Resumen eigen-espacio homólogo corregido:\n")
cat("  NaN      :", num_nan, "\n")
cat("  Inf      :", num_inf, "\n")
cat("  Imaginarios distintos de 0:", tiene_parte_imaginaria, "\n")
cat("  Mínimo   :", min_val, "\n")
cat("  Máximo   :", max_val, "\n")

# Criterio de seguridad:
# - si hay NaN/Inf/imaginarios -> detener
# - si los valores son absurdamente grandes (> 2, por ejemplo), detener
if (num_nan > 0 || num_inf > 0 || tiene_parte_imaginaria) {
  stop("Los eigenvectores corregidos contienen NaN, Inf o parte imaginaria. Revisa antes de continuar.")
}

if (max_val > 2) {
  stop("Los eigenvectores corregidos tienen valores > 2, lo cual es sospechoso para loadings de PCA. Revisa la transformación.")
}

## =========================
## 9) Guardar eigen-espacio corregido
## =========================
saveRDS(
  loadings_homolog_corregidos,
  file = file.path(OUT_DIR, "Ovary_eigenspace_homolog_loadings_corrected.rds")
)
