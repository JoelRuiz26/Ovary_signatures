# =============================================================================
# Master Regulator Analysis con paquete corto
# Usando los datos de ejemplo del propio paquete (inmat y centroids)
# =============================================================================

library(corto)

# -----------------------------------------------------------------------------
# 1. Cargar datos de ejemplo del paquete
# -----------------------------------------------------------------------------
load(system.file("extdata", "inmat.rda", package = "corto", mustWork = TRUE))
load(system.file("extdata", "centroids.rda", package = "corto", mustWork = TRUE))

# inmat: matriz de expresión (7000 genes x 87 muestras)
# centroids: vector con 294 factores de transcripción

# -----------------------------------------------------------------------------
# 2. Simular dos condiciones: "relapse" y "control"
#    Dividimos las 87 muestras en dos grupos (mitad y mitad)
# -----------------------------------------------------------------------------
n_total <- ncol(inmat)
n_control <- floor(n_total / 2)          # 43 muestras
n_relapse <- n_total - n_control         # 44 muestras

counts_control <- inmat[, 1:n_control]
counts_relapse <- inmat[, (n_control+1):n_total]

# Filtro opcional: asegurar que los centroids existan en ambas matrices
tfs <- centroids   # usamos los centroids del paquete como lista de TFs

# -----------------------------------------------------------------------------
# 3. Construir regulones para cada condición (opcional, pero como en tu script)
# -----------------------------------------------------------------------------
centroids_control <- tfs[tfs %in% rownames(counts_control)]
regulon_control <- corto(counts_control, 
                         centroids = centroids_control, 
                         nbootstraps = 100,   # en ejemplo rápido podrías bajar a 10
                         p = 1e-3, 
                         nthreads = 6,
                         verbose = TRUE)

centroids_relapse <- tfs[tfs %in% rownames(counts_relapse)]
regulon_relapse <- corto(counts_relapse, 
                         centroids = centroids_relapse, 
                         nbootstraps = 100, 
                         p = 1e-3, 
                         nthreads = 6,
                         verbose = TRUE)



# -----------------------------------------------------------------------------
# 4. Análisis de expresión diferencial con DESeq2
# -----------------------------------------------------------------------------

# Asegurar que las matrices de conteo sean enteros (DESeq2 requiere integers)
counts_relapse <- round(counts_relapse)
counts_control <- round(counts_control)

# 4.1 Unir las dos matrices en una sola (columnas: controles + recaídas)
counts_combined <- cbind(counts_control, counts_relapse)

# 4.2 Crear el data.frame con la información de las muestras (colData)
colData <- data.frame(
  condition = factor(c(rep("control", ncol(counts_control)),
                       rep("relapse", ncol(counts_relapse)))),
  row.names = colnames(counts_combined)
)

# 4.3 Crear el objeto DESeqDataSet
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_combined,
                              colData = colData,
                              design = ~ condition)

# 4.4 Pre‑filtrado opcional (elimina genes con muy pocas lecturas)
keep <- rowSums(counts(dds) >= 10) >= 3   # al menos 3 muestras con cuenta ≥ 10
dds <- dds[keep, ]

# 4.5 Ejecutar DESeq2 (estimación de dispersión, GLM, test de Wald)
dds <- DESeq(dds)

# 4.6 Extraer resultados: comparación relapse vs control
res <- results(dds, contrast = c("condition", "relapse", "control"))

# 4.7 (Opcional) Ordenar por p‑valor y visualizar primeros resultados
resOrdered <- res[order(res$pvalue), ]
head(resOrdered)

# -----------------------------------------------------------------------------
# 5. Crear firma de expresión diferencial (estadístico t)
# -----------------------------------------------------------------------------
t_stat <- res$stat          # ya es el estadístico de Wald
names(t_stat) <- rownames(res)

# Eliminar NA (si los hay)
t_stat <- t_stat[!is.na(t_stat) & is.finite(t_stat)]

# Filtrar por genes presentes en tu regulón
genes_regulon <- unique(unlist(sapply(regulon_relapse, function(x) names(x$tfmode))))
t_stat <- t_stat[names(t_stat) %in% genes_regulon]


# -----------------------------------------------------------------------------
# 6. Master Regulator Analysis usando la firma de DESeq2
#    (Requiere la función modificada mra_mod que admite vectores)
# -----------------------------------------------------------------------------
source("~/Ovary_signatures/7_Survival_analysis/mra_mod.R")
mrs <- mra_mod(expmat1 = t_stat,          # firma precomputada
               expmat2 = NULL,
               regulon = regulon_relapse,
               minsize = 20,
               nperm = 1000,
               nthreads = 6,
               verbose = TRUE)




# -----------------------------------------------------------------------------
# 5. Gráfico con mraplot (guarda en PDF)
# -----------------------------------------------------------------------------
pdf(file = "cortoMRS_ejemplo.pdf", width = 15, height = 18)
mraplot(mrs, mrs = 10)
dev.off()

# -----------------------------------------------------------------------------
# 6. Tabla de resultados con NES, p-value y FDR
# -----------------------------------------------------------------------------
alls <- data.frame(TF = names(mrs$nes),
                   NES = mrs$nes,
                   pvalue = mrs$pvalue)
alls$FDR <- p.adjust(alls$pvalue, method = "fdr", n = length(alls$pvalue))

# Top reguladores con p < 0.05
tops <- alls[alls$pvalue < 0.05, ]
tops <- tops[order(tops$pvalue), ]

# Mostrar resultados
cat("\n=== Top 10 Master Regulators (p < 0.05) ===\n")
print(head(tops, 10))

# Guardar tabla completa
write.csv(alls, "master_regulators_completos.csv", row.names = FALSE)
write.csv(tops, "master_regulators_significativos.csv", row.names = FALSE)

cat("\nResultados guardados en 'master_regulators_completos.csv' y 'master_regulators_significativos.csv'\n")