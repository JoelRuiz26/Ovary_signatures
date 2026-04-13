
set JAVA_TOOL_OPTIONS "-XX:ActiveProcessorCount=30"
java -Xmx50G -jar /STORAGE/genut/software/ARACNe-AP/dist/aracne.jar -e ov_tcga_vst.tsv -o cancer_ovary_network_300bt_p1e-8 --tfs centroids/Regulators_All_TF_TFC_miRNA.txt --pvalue 1E-8 --seed 1 --calculateThreshold

for i in (seq 1 300)
  echo "Run $i/300"
  java -Xmx200G -jar /STORAGE/genut/software/ARACNe-AP/dist/aracne.jar -e ov_tcga_vst.tsv -o cancer_ovary_network_300bt_p1e-8 --tfs centroids/Regulators_All_TF_TFC_miRNA.txt --pvalue 1E-8 --threads 50 --seed $i
end

java -Xmx50G -jar /STORAGE/genut/software/ARACNe-AP/dist/aracne.jar -o cancer_ovary_network_300bt_p1e-8 --consolidate

mv cancer_ovary_network_300bt_p1e-8/network.txt ./cancer_ovary_network_300bt_p1e-8.txt

