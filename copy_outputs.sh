#!/bin/bash
set -e

echo "=== Sincronizando Output 1 (RawOvary) ==="
rsync -av /home/jruiz/Ovary_signatures/1_DGE_signature_ovary/1_1_Output_rds/ \
          /STORAGE/csbig/jruiz/Ovary_data/1_Output_RawOvary/

echo "=== Sincronizando Output 2 (AutoEncoder - adjacent) ==="
# Crear carpeta si no existe
mkdir -p /STORAGE/csbig/jruiz/Ovary_data/2_Output_AutoEncoder/adjacent/
rsync -av /home/jruiz/Ovary_signatures/2_DGE_signature_autoEncoder/2_1_Output_adjacent_rds/ \
          /STORAGE/csbig/jruiz/Ovary_data/2_Output_AutoEncoder/adjacent/

echo "=== Sincronizando Output 2 (AutoEncoder - gtex) ==="
# Crear carpeta si no existe
mkdir -p /STORAGE/csbig/jruiz/Ovary_data/2_Output_AutoEncoder/gtex/
rsync -av /home/jruiz/Ovary_signatures/2_DGE_signature_autoEncoder/2_1_Output_gtex_rds/ \
          /STORAGE/csbig/jruiz/Ovary_data/2_Output_AutoEncoder/gtex/

echo "=== Sincronizando Output 3 (FixedOvary) ==="
rsync -av /home/jruiz/Ovary_signatures/3_DGE_signature_correctedOvary/ \
          /STORAGE/csbig/jruiz/Ovary_data/3_Output_FIxedOvary/

echo "=== Sincronización completa ==="
