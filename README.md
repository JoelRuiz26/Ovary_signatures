# Ovary_signatures — Summary

## 0. 0_RawCounts_Ovary
Matriz cruda de conteos (TCGA ovary + GTEx ovary + autoencoder). Punto de partida para todas las firmas posteriores.

## 1. 1_DGE_signature_ovary
Expresión diferencial entre TCGA ovarian cancer vs GTEx ovary homologous.
Incluye:
- Script principal (1_0_DiseaseSignature_3m_all.R)
- Resultados DE (DESeq2, EdgeR, limma)
- Meta-pValues (Stouffer)
- Metadata y tablas finales (1_1_Output_rds)

## 2. 2_DGE_signature_autoEncoder
Expresión diferencial entre TCGA ovarian cancer y controles seleccionados por autoencoder.

### adjacent/
Controles autoencoder provenientes de GTEx + TCGA adjacent + TARGET.

### gtex/
Controles autoencoder solo de GTEx.

Ambos incluyen:
- Scripts (2_0_DGE_3Methods_adjacent.R / 2_0_DGE_3Methods_gtex.R)
- Resultados DE (3 métodos)
- Meta-pValues
- Volcano plot de la metasignature

## 3. 3_DGE_signature_correctedOvary
Comparación entre TCGA ovarian cancer y GTEx ovary corregido (raw counts filtrados + QC + corrección de sesgos).
Incluye:
- Scripts de construcción de matrices corregidas
- Metadata filtrada
- PCA, t-SNE y UMAP
- Raw counts ajustados para autoencoder y homologous (subcarpetas adjacent/ y gtex/)

