# =========================================
# Load Required Libraries
# =========================================
library(octad)
library(octad.db)
library(dplyr)
library(vroom)
library(tibble)
library(ggplot2)
library(harmonicmeanp)  


setwd("~/Ovary_signatures/2_DGE_signature_autoEncoder/")

# ==
### ————— After running the 3 methods, load the three results and generate the meta-signature:
### 1. Load results
edge   <- readRDS("2_1_Output_adjacent_rds/DE_OVARY_ALL_EdgeR.rds")
deseq2 <- readRDS("2_1_Output_adjacent_rds/DE_OVARY_ALL_DESeq2.rds")
limma  <- readRDS("2_1_Output_adjacent_rds/DE_OVARY_ALL_limma.rds")

### 2. Merge tables (logFC + p-values) 
# https://www.pnas.org/doi/full/10.1073/pnas.1814092116
# Para combinar la evidencia de los tres métodos de expresión diferencial
# (edgeR, DESeq2 y limma), utilizamos el método inverse normal / Stouffer.
# Primero convertimos los p-values de cada método en Z-scores con signo,
# usando el signo del logFC. Después combinamos los Z por gen:
#   Z_meta = sum(Z_i) / sqrt(k) (k=nMethods)
# y obtenemos un meta p-value bilateral:
#   p_meta = 2 * pnorm(-abs(Z_meta)).
###pnorm() es la función de distribución acumulada de la normal estándar N(0,1)
###es decir, probabilidad de que una variable normal estándar sea ≤ z
# Finalmente, corregimos por múltiples pruebas con Benjamini–Hochberg.

meta_tmp <- edge %>%
  dplyr::select(
    Symbol,
    log2FC_edgeR = log2FoldChange,
    p_edgeR      = pvalue,
    padj_edgeR   = padj
  ) %>%
  dplyr::inner_join(
    deseq2 %>% dplyr::select(
      Symbol,
      log2FC_DESeq2 = log2FoldChange,
      p_DESeq2      = pvalue,
      padj_DESeq2   = padj
    ),
    by = "Symbol"
  ) %>%
  dplyr::inner_join(
    limma %>% dplyr::select(
      Symbol,
      log2FC_limma = log2FoldChange,
      p_limma      = pvalue,
      padj_limma   = padj
    ),
    by = "Symbol"
  )

### 3. Construir matrices de p-values y logFC (genes x métodos)
pmat <- as.matrix(meta_tmp[, c("p_edgeR", "p_DESeq2", "p_limma")])
lfc_mat <- as.matrix(meta_tmp[, c("log2FC_edgeR", "log2FC_DESeq2", "log2FC_limma")])

# Evitar problemas numéricos con p extremadamente pequeños o cercanos a 1
eps <- .Machine$double.xmin
pmat[pmat < eps]       <- eps
pmat[pmat > 1 - eps]   <- 1 - eps

### 4. Calcular Z-scores con signo para cada método
# p de entrada se asumen bilaterales; construimos Z firmado:
# Z_i = sign(logFC_i) * qnorm(p_i/2, lower.tail = FALSE)
z_mat <- sign(lfc_mat) * qnorm(pmat / 2, lower.tail = FALSE)

### 5. Z combinado por gen (Stouffer) y meta p-value bilateral
k_methods <- ncol(z_mat)

z_meta <- rowSums(z_mat, na.rm = TRUE) / sqrt(k_methods)
p_meta <- 2 * pnorm(-abs(z_meta))  # bilateral
p_meta[p_meta < eps] <- eps


### 6. Construir meta_all con Stouffer
meta_all <- meta_tmp %>%
  dplyr::mutate(
    z_meta       = z_meta,
    pvalue_meta  = p_meta,
    padj_meta    = p.adjust(pvalue_meta, method = "BH"),
    log2FC_mean  = (log2FC_edgeR + log2FC_DESeq2 + log2FC_limma) / 3,
    regulation   = dplyr::if_else(log2FC_mean > 0, "up", "down"),
    # cuántos métodos individuales son significativos (FDR < 0.01)
    n_sig_methods = (padj_edgeR  < 0.01) +
      (padj_DESeq2 < 0.01) +
      (padj_limma  < 0.01)
  ) %>%
  dplyr::select(
    Symbol,
    log2FoldChange = log2FC_mean,
    regulation,
    z_meta,
    pvalue_meta,
    padj = padj_meta,
    n_sig_methods,
    log2FC_edgeR,  p_edgeR,  padj_edgeR,
    log2FC_DESeq2, p_DESeq2, padj_DESeq2,
    log2FC_limma,  p_limma,  padj_limma
  )
meta_all$padj[meta_all$padj < eps] <- eps


saveRDS(
  meta_all,
  "2_1_Output_adjacent_rds/2_2_Meta_DE_OVARY_ALL_Stouffer.rds"
)

### 7. Meta-significant & consistent genes (consenso + efecto)
meta_sig <- meta_all %>%
  dplyr::filter(
    # consistencia de dirección entre métodos
    sign(log2FC_edgeR)   == sign(log2FC_DESeq2),
    sign(log2FC_edgeR)   == sign(log2FC_limma),
    # significancia según p combinada (Stouffer)
    padj <= 0.01,
    # al menos 2 métodos individuales significativos (puedes cambiar a ==3 si quieres ultra estricto)
    n_sig_methods >= 2,
    # tamaño de efecto mínimo
    abs(log2FoldChange) >= 1
  )

saveRDS(
  meta_sig,
  "2_1_Output_adjacent_rds/2_2_Meta_DE_OVARY_significant_Stouffer.rds"
)