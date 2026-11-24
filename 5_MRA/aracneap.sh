#!/usr/bin/env bash
set -euo pipefail

# ARACNe-AP network of TCGA-OV tumors: 300 bootstraps, MI p-value 1e-8.
# Run from the folder that holds ov_tcga_vst.tsv and
# modulators_all_tf_tfc_mirna.txt (both written by 5_00_expression_regulators.R).
ARACNE_JAR="${ARACNE_JAR:-aracne.jar}"   # https://github.com/califano-lab/ARACNe-AP
THREADS="${ARACNE_THREADS:-50}"

export JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=30"
java -Xmx50G -jar "$ARACNE_JAR" -e ov_tcga_vst.tsv -o cancer_ovary_network_300bt_p1e-8 --tfs modulators_all_tf_tfc_mirna.txt --pvalue 1E-8 --seed 1 --calculateThreshold

for i in $(seq 1 300); do
  echo "Run $i/300"
  java -Xmx200G -jar "$ARACNE_JAR" -e ov_tcga_vst.tsv -o cancer_ovary_network_300bt_p1e-8 --tfs modulators_all_tf_tfc_mirna.txt --pvalue 1E-8 --threads "$THREADS" --seed "$i"
done

java -Xmx50G -jar "$ARACNE_JAR" -o cancer_ovary_network_300bt_p1e-8 --consolidate

mv cancer_ovary_network_300bt_p1e-8/network.txt ./cancer_ovary_network_300bt_p1e-8.txt
