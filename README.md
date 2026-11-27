# Ovary Signatures

R pipeline supporting the manuscript:

> Ruiz-Hernández et al. *Integrative transcriptomic and network-guided functional dependency analysis reveals essential transcriptional drivers in ovarian cancer* (under revision).

The workflow has three layers: (1) a consensus differential-expression signature of ovarian serous cystadenocarcinoma from three sources (GEO microarrays, TCGA vs GTEx ovary, TCGA vs autoencoder-selected controls); (2) Master Regulator Analysis (ARACNe-AP + msVIPER) over transcription factors and cofactors; (3) prioritization of regulators with CRISPR dependency data from DepMap.

![Graphical abstract](assets/graphical_abstract.png)

## Results in brief

- Consensus signature: 698 genes significant in the three sources (382 up, 316 down; direction by agreement of at least two sources).
- Master regulators: 2,436 regulons tested, 1,469 significant (BH-adjusted empirical p < 0.05); 227 above the elbow of the |NES| curve.
- DepMap 24Q4, 58 ovarian epithelial tumor lines: 64 consensus genes with median gene effect ≤ −1; 69 regulators enriched in them (Fisher test against the consensus genes with dependency data); 21 candidate regulators that are themselves dependencies (median gene effect ≤ −0.5) and expressed in the same lines (median TPM ≥ 1): ACTL6A, AURKB, BIRC5, CCNA2, CDK1, CDK2, CENPA, CHAF1B, CHEK1, DNMT1, GMNN, H2AZ1, ILF2, MIS18BP1, PLK1, PTTG1, SRSF2, TAF5, TIMELESS, TTF2, WDHD1.

## Layout and paths

The code lives in this repository; data and results live under a **data root** with the same folder names (`0_DGE_GTEx/`, `1_DGE_AE/`, …). Every script reads the data root from the environment variable `OVARY_SIGNATURES` (default `~/Ovary_signatures`, i.e. the repository itself used as data root):

```bash
export OVARY_SIGNATURES=/path/to/data_root
Rscript 5_MRA/5_0_ovary_mra.R
```

Other variables: `OCTAD_H5` (path to `octad.counts.and.tpm.h5`), `DEPMAP_DIR` (default `$OVARY_SIGNATURES/data/depmap_24Q4`), `ARACNE_JAR`, and the `MRA_*` settings documented in the header of `5_MRA/5_0_ovary_mra.R`.

## External inputs

| Input | Used by | Where to get it |
|---|---|---|
| `octad.counts.and.tpm.h5` | `0_1`, `1_1`, `3_0` | OCTAD expression file (see the `octad` package documentation) |
| GEO series (six microarray datasets) | `2_DEG_GEO/GEO_multi_RRA.R` | downloaded by the script with `GEOquery` |
| `TcgaTargetGtex_gene_expected_count.gz`, `TcgaTargetGTEX_phenotype.txt.gz` | `5_00` | UCSC Xena Toil hub, `https://toil.xenahubs.net/download/` → `data/xena/` |
| `Homo_sapiens_TF.txt`, `Homo_sapiens_Cof.txt` | `5_00` | AnimalTFDB → `data/regulators/` |
| `DatabaseExtract_v_1.01.txt` | `5_00` | Lambert et al. 2018, `https://humantfs.ccbr.utoronto.ca/download.php` → `data/regulators/` |
| `CRISPRGeneEffect.csv`, `Model.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | `6_*`, `7_2` | DepMap Public 24Q4, `https://depmap.org/portal/data_page/` → `data/depmap_24Q4/` |
| ARACNe-AP (`aracne.jar`) | `5_MRA/aracneap.sh` | `https://github.com/califano-lab/ARACNe-AP` |

## Order of execution

| Step | Script | Main outputs |
|---|---|---|
| 1 | `0_DGE_GTEx/0_1_Deseq_octad_GTEX.R` | `DE_full_OVARY_DESeq2_GTEx.rds` |
| 2 | `1_DGE_AE/1_1_Deseq_octad_AE.R` | `DE_full_OVARY_DESeq2_AE.rds` |
| 3 | `2_DEG_GEO/GEO_multi_RRA.R` | `GEO_ovarian_cancer_RRA_results.rds`, `GEO_ovarian_cancer_RRA_1row.rds` |
| 4 | `3_Consensus_DGE_analysis/3_0_Counts_norm.R` | normalized counts, metadata, gene annotation (`3_0_Counts_metadata/`) |
| 5 | `3_Consensus_DGE_analysis/3_1_DGE_list_shared.R` | consensus tables (`3_1_Core_signature/3_0_2_allgenes.tsv`), Venn diagram |
| 6 | `3_2_Heatmap.R`, `3_2_Volcano_plots_core_signature.R`, `3_3_Fano_poison.R` | heatmaps, volcano plots, dispersion analysis |
| 7 | `4_ORA_clusterProfiler/4_ORA_core_bar.R` | GO over-representation of the consensus |
| 8 | `5_MRA/5_00_expression_regulators.R` | `ov_tcga_vst.tsv`, `modulators_all_tf_tfc_mirna.txt`, `onlyTFs.txt` |
| 9 | `5_MRA/aracneap.sh` | `cancer_ovary_network_300bt_p1e-8.txt` |
| 10 | `5_MRA/5_0_ovary_mra.R` | Figures S1–S2, `signature_core.tsv`, regulon object, `core_mra.tsv` |
| 11 | `5_MRA/5_1_Elbow_TF.R` | top regulators by \|NES\| (elbow) |
| 12 | `6_Depmap_ovary/6_1_Heatmap_coregenes.R` → `6_2_MR_Enriched.R` → `6_3_Heatmap_TF.R` | dependency of consensus genes, enriched regulators, candidate drivers |
| 13 | `5_MRA/5_2_Circos_regulon_overlap.R`, `5_MRA/5_3_Circos_driver_hallmark.R` | circos plots of the candidates (need the raw network) |
| 14 | `7_Figures/7_1_…`, `7_2_…`, `7_3_…` | assembled manuscript figures (`7_3` reads the caches written by `7_2`) |
| 15 | `8_Survival_analysis/8_0_Clinical_metadata.R` → `8_1_Cox_drivers_survival.R` | TCGA-OV clinical table, Cox models (`8_1` reads the caches written by `7_2`) |

In step 3 the GEO cases are high-grade serous carcinomas and tumor epithelium only, and the controls are ovarian surface epithelium (low-grade and borderline tumors and microdissected stroma are left out). The probe-to-symbol mapping is read from `2_DEG_GEO/annotation/` (frozen from g:Profiler; version and date in the file header). Step 8 reads the gene annotation from the dictionaries frozen in `5_MRA/annotation/` (Ensembl release 115; each file records release and date in its header), so it does not depend on the Ensembl servers; delete a dictionary to rebuild it with biomaRt. With these inputs `5_00` writes a matrix that is byte-identical to the one used for the network, and an ARACNe-AP bootstrap run with the same seed reproduces the corresponding bootstrap file of the original run. `aracneap.sh` must be run from `$OVARY_SIGNATURES/5_MRA`.

Steps 9 and 10 are the expensive ones: the network takes about 6.5 h with 50 threads, and 10,000 permutations of msVIPER take a few minutes per signature with 55 workers. The permutation null of `5_0` uses one seed per permutation, so the result does not depend on the number of workers.

## Running with Snakemake

The `Snakefile` only encodes the order of execution: one rule per script (steps 1 to 15 and `9_Revision_analyses/`) with its main inputs and outputs. The scripts run unchanged and read the same environment variables.

```bash
export OVARY_SIGNATURES=/path/to/data_root
snakemake -n            # list what is missing or older than its inputs (file timestamps)
snakemake --cores 1     # run it, one script at a time; the targets `pipeline` and `revision` run one part only
```

## Revision analyses

The scripts in `9_Revision_analyses/` are the additional analyses of the revision. They read the outputs of the pipeline from the data root and write to `$OVARY_SIGNATURES/9_Revision_analyses/<subfolder>/`. The candidate regulators are always read from `6_Depmap_ovary/6_3_2_candidate_regulators.tsv`, and the elbow from `5_MRA/5_0_1_Top_elbow_regulones.rds`.

| Script | Main outputs | Question |
|---|---|---|
| `9_1_1_Fano_tests.R` | `consensus/fano_tests.tsv`, `fano_paired_cohorts.tsv` | Is the Fano factor of the consensus genes different from that of the other differentially expressed genes, also after removing the mean–variance trend? |
| `9_1_2_ORA_by_direction.R` | `consensus/ORA_GO_consensus_by_direction.tsv` | Which GO terms are over-represented among the up- and the down-regulated consensus genes separately? |
| `9_2_1_Prioritization_sensitivity.R` | `prioritization/prioritization_sensitivity.tsv`, `null_random_HD_consensus*.tsv`, `fisher_by_regulator.tsv`, `HD_targets_by_regulator.tsv` | How do the candidates change when each threshold of the prioritization is varied, and what do random sets of dependent genes give? |
| `9_2_2_Lineage_selectivity.R` | `prioritization/lineage_selectivity_candidates.tsv`, `ovarian_cell_lines.tsv`, `Figure_S_lineage` | Are the dependencies on the candidates selective for ovarian lines or common to all lineages? |
| `9_3_1_Export_regulons.R` → `9_3_2_Perturbseq_regulon_response.py` → `9_3_3_Figure_perturbseq.R` | `perturbseq/perturbseq_regulon_tests.tsv`, `Figure_S_perturbseq` | Do the positive targets of each regulon go down when the regulator is knocked down (Perturb-seq)? |
| `9_4_1_FT_series_metadata.R` → `9_4_2_FT_limma.R` → `9_4_3_FT_consensus_concordance.R` → `9_4_4_FT_mra.R` → `9_4_5_FT_activated_vs_repressed.R` → `9_4_6_FT_figure_NES.R` | `fallopian_tube/de_*.tsv`, `consensus_concordance.tsv`, `ft_mra_vs_main.tsv`, `ft_mra_candidates.tsv`, `elbow_activated_vs_repressed.tsv`, `Figure_S_NES_main_vs_fallopian_tube` | Do the consensus signature and the master regulators hold when tumors are compared with fallopian tube epithelium instead of ovary? |
| `9_5_1_Cox_meta_independent.R` → `9_5_2_Forest_independent.R` | `survival_independent/meta_by_gene.tsv`, `cox_by_dataset.tsv`, `group_tests.tsv`, `comparison_TCGA.tsv`, `forest_meta_independent` | Is the expression of the candidates associated with overall survival in cohorts independent of TCGA? |
| `9_6_1_Octad_reference_metric.R` → `9_6_2_Octad_reference_figure.R` | `octad_metric/metric_by_sample.tsv`, `group_summary.tsv`, `composition_AE_controls.tsv`, `octad_reference_metric` | Which similarity metric selects the autoencoder controls, and where do the GTEx ovary controls fall on it? |
| `9_7_1_HPA_download_xml.sh` → `9_7_2_HPA_extract.py` → `9_7_3_HPA_pattern_rule.py` | `hpa/hpa_v25_ovary_by_antibody.tsv`, `hpa_v25_ovary_patients.tsv`, `pattern_rule_by_antibody.tsv` | Is the protein absent or low in normal ovarian stroma and detected in most ovarian tumors (immunohistochemistry)? |

Additional external inputs:

| Input | Used by | Where to get it |
|---|---|---|
| `K562_essential_normalized_bulk_01.h5ad`, `K562_gwps_normalized_bulk_01.h5ad`, `rpe1_normalized_bulk_01.h5ad` | `9_3_2` | Replogle et al. 2022, processed Perturb-seq datasets, `https://doi.org/10.25452/figshare.plus.20029387` → `data/replogle2022/` (or `REPLOGLE_DIR`) |
| GEO series GSE10971 and GSE69428 (GPL570) | `9_4_1`, `9_4_2` | downloaded by `9_4_1` with `GEOquery` → `9_Revision_analyses/fallopian_tube/geo_data/` |
| `curatedOvarianData` | `9_5_1` | Bioconductor package (Ganzfried et al. 2013) |
| Human Protein Atlas v25, one XML per gene | `9_7_2` | downloaded by `9_7_1` from `https://v25.proteinatlas.org/` for the genes of `9_7_0_HPA_genes.tsv` → `9_Revision_analyses/hpa/xml/` |

`9_4_2` uses the probe annotation frozen in `2_DEG_GEO/annotation/`. `9_6_1` needs the `octad` package, as steps 1, 2 and 4. With 10 workers `9_4_4` (1,000 permutations of msVIPER per signature) takes about 17 minutes and `9_5_1` (10,000 random gene sets and 10,000 permutations) about 3; the random draws of `9_5_1` depend on the number of workers, which is fixed in the script.

## Dependencies

R ≥ 4.4 with: `tidyverse`, `vroom`, `janitor`, `DESeq2`, `limma`, `octad`, `octad.db`, `GEOquery`, `RobustRankAggreg`, `gprofiler2`, `biomaRt`, `clusterProfiler`, `enrichplot`, `org.Hs.eg.db`, `msigdbr`, `viper`, `ComplexHeatmap`, `pheatmap`, `dynamicTreeCut`, `circlize`, `ggvenn`, `ggrepel`, `ggh4x`, `patchwork`, `png`, `effsize`, `TCGAbiolinks`, `survival`, `data.table`, `curatedOvarianData`, `ragg`. Java for ARACNe-AP. Python ≥ 3.10 with `anndata`, `pandas`, `numpy` and `scipy` ≥ 1.11 for the Perturb-seq script (the Human Protein Atlas scripts use the standard library and `curl`). Snakemake is optional (tested with 9.14).

## License

GPL-3

## Contact

Open a GitHub issue, or contact the corresponding authors listed in the manuscript.
