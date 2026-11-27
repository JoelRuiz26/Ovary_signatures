# Order of execution of the pipeline and of the revision analyses.
#
# Thin wrapper: one rule per script, with its main inputs and outputs, and nothing
# else. The scripts are run as they are and take every setting from the environment
# (OVARY_SIGNATURES, DEPMAP_DIR, OCTAD_H5, ARACNE_JAR, MRA_*; see README.md).
#
#   export OVARY_SIGNATURES=/path/to/data_root
#   snakemake -n                # what would run
#   snakemake --cores 1         # run whatever is missing or out of date
#
# --cores only limits how many scripts run at once; the threads used inside a script
# are set by its own variables (MRA_CORES, ARACNE_THREADS).

import os

# Data folders are kept with a trailing slash, so that paths are built as FOLDER + "file"
def folder(path):
    return os.path.join(os.path.expanduser(path), "")

REPO     = workflow.basedir                              # scripts of the pipeline
RSRC     = os.path.join(REPO, "9_Revision_analyses")     # scripts of the revision analyses
ROOT     = folder(os.environ.get("OVARY_SIGNATURES", "~/Ovary_signatures"))
DEPMAP   = folder(os.environ.get("DEPMAP_DIR", ROOT + "data/depmap_24Q4"))
REPLOGLE = folder(os.environ.get("REPLOGLE_DIR", ROOT + "data/replogle2022"))
# OCTAD count file: declared as an input only when the variable is set
OCTAD_H5 = [os.environ["OCTAD_H5"]] if "OCTAD_H5" in os.environ else []

DGE  = ROOT + "3_Consensus_DGE_analysis/"
CNT  = DGE + "3_0_Counts_metadata/"
MRA  = ROOT + "5_MRA/"
DEP  = ROOT + "6_Depmap_ovary/"
FIG  = ROOT + "7_Figures/"
SURV = ROOT + "8_Survival_analysis/"
REV  = ROOT + "9_Revision_analyses/"

DE_GTEX    = ROOT + "0_DGE_GTEx/DE_full_OVARY_DESeq2_GTEx.rds"
DE_AE      = ROOT + "1_DGE_AE/DE_full_OVARY_DESeq2_AE.rds"
GEO_RRA    = ROOT + "2_DEG_GEO/GEO_ovarian_cancer_RRA_results.rds"
GEO_1ROW   = ROOT + "2_DEG_GEO/GEO_ovarian_cancer_RRA_1row.rds"
ALLGENES   = DGE + "3_1_Core_signature/3_0_2_allgenes.tsv"
CONCORDANT = DGE + "3_1_Core_signature/3_0_1_genes_concordant.tsv"
FANO       = expand(DGE + "3_3_Fano_Factor/3_3_3_fano_{c}_raw.rds", c=["GTEX", "AE"])
ORA        = ROOT + "4_ORA_clusterProfiler/ORA_GO_consensus3of3_ALL_BP_MF_CC_SIMPLIFIED_0.3.tsv"
VST        = MRA + "ov_tcga_vst.tsv"
REGULATORS = MRA + "modulators_all_tf_tfc_mirna.txt"
ONLY_TFS   = MRA + "onlyTFs.txt"
NETWORK    = MRA + "cancer_ovary_network_300bt_p1e-8.txt"
CORE_MRA   = MRA + "core_mra.tsv"
SIGNATURE  = MRA + "signature_core.tsv"
ELBOW_TAB  = MRA + "5_0_1_Top_elbow_regulones.rds"
ELBOW_REG  = MRA + "5_0_Top_elbow_regulones.rds"
DEP_CORE   = DEP + "6_1_0_Dependency_DEG_core.tsv"
CANDIDATES = DEP + "6_3_2_candidate_regulators.tsv"
MRA_CACHE  = FIG + "7_2_mra_enrichment_24Q4_646_cache.rds"
COX_TCGA   = SURV + "8_1_0_cox_drivers_survival.tsv"
CRISPR     = [DEPMAP + "Model.csv", DEPMAP + "CRISPRGeneEffect.csv"]
DEPMAP_TPM = DEPMAP + "OmicsExpressionProteinCodingGenesTPMLogp1.csv"

# The regulon object (5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds) and the dependency cache of 7_2
# (7_Figures/7_2_depmap_ovary_24Q4_646_cache.rds) are rewritten by their scripts only when they are
# older than the files they are built from, so their timestamps say nothing about the other inputs
# of the rule. They are not declared as outputs: the file that the same script always writes after
# them stands for them wherever a rule reads them.
REGULON    = CORE_MRA    # 5_0_ovary_mra.R writes core_mra.tsv after the regulon object
DEP_CACHE  = MRA_CACHE   # 7_2 writes the enrichment cache after the dependency cache

HPA_GENES  = [l.split()[1] for l in open(os.path.join(RSRC, "9_7_0_HPA_genes.tsv")) if l.strip()]


# ---------------------------------------------------------------------------
# Main pipeline (steps 1 to 15 of the README)
# ---------------------------------------------------------------------------

# Step 1
rule dge_gtex:
    input: OCTAD_H5
    output: DE_GTEX
    shell: "Rscript {REPO}/0_DGE_GTEx/0_1_Deseq_octad_GTEX.R"

# Step 2
rule dge_ae:
    input: OCTAD_H5
    output: DE_AE
    shell: "Rscript {REPO}/1_DGE_AE/1_1_Deseq_octad_AE.R"

# Step 3 (downloads the GEO series; its last section reads the two DESeq2 tables)
rule geo_rra:
    input: DE_GTEX, DE_AE
    output: GEO_RRA, GEO_1ROW, ROOT + "2_DEG_GEO/DEGs_alldsets_GEO.rds"
    shell: "Rscript {REPO}/2_DEG_GEO/GEO_multi_RRA.R"

# Step 4
rule counts_norm:
    input: OCTAD_H5
    output:
        expand(CNT + "{f}.rds", f=["3_3_1_expr_raw_Tumor_Ctlhomolog", "3_3_1_expr_raw_Tumor_CtlAuto",
                                    "3_3_1_metadata_all", "3_3_2_dds_tumor_Ctlhomolog", "3_3_2_dds_tumor_CtlAuto",
                                    "3_3_2_norm_counts_tumor_GTEX", "3_3_2_norm_count_tumor_AE",
                                    "3_3_2_norm_counts_tumor_ONLY", "3_3_2_gene_annotation_UNION_universe"])
    shell: "Rscript {REPO}/3_Consensus_DGE_analysis/3_0_Counts_norm.R"

# Step 5
rule consensus:
    input: DE_GTEX, DE_AE, GEO_RRA, CNT + "3_3_2_gene_annotation_UNION_universe.rds"
    output: ALLGENES, CONCORDANT, expand(DGE + "3_1_Core_signature/3_0_1_venn_TCGA_GEO.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/3_Consensus_DGE_analysis/3_1_DGE_list_shared.R"

# Step 6
rule consensus_heatmaps:
    input: CNT + "3_3_2_dds_tumor_Ctlhomolog.rds", CNT + "3_3_2_dds_tumor_CtlAuto.rds", ALLGENES
    output:
        expand(DGE + "3_2_DGE_volcano_plots/3_2_2_HM_core_{c}_all.png", c=["GTEX", "AE"]),
        expand(DGE + "3_2_DGE_volcano_plots/3_2_2_HM_clusters_{c}_all.tsv", c=["GTEX", "AE"])
    shell: "Rscript {REPO}/3_Consensus_DGE_analysis/3_2_Heatmap.R"

rule consensus_volcano:
    input: DE_GTEX, DE_AE, CONCORDANT
    output: expand(DGE + "3_2_DGE_volcano_plots/3_2_1_volcano_GTEx_AE_only_CORE.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/3_Consensus_DGE_analysis/3_2_Volcano_plots_core_signature.R"

rule fano:
    input: CNT + "3_3_2_norm_counts_tumor_GTEX.rds", CNT + "3_3_2_norm_count_tumor_AE.rds", ALLGENES
    output: FANO, expand(DGE + "3_3_Fano_Factor/3_3_4_fano_subset_lines.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/3_Consensus_DGE_analysis/3_3_Fano_poison.R"

# Step 7
rule ora:
    input: DE_GTEX, DE_AE, GEO_RRA
    output: ORA, expand(ROOT + "4_ORA_clusterProfiler/ORA_GO_topTerms_dotplot_COMPACT_ALL.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/4_ORA_clusterProfiler/4_ORA_core_bar.R"

# Step 8
rule expression_regulators:
    input:
        ROOT + "data/xena/TcgaTargetGtex_gene_expected_count.gz",
        ROOT + "data/xena/TcgaTargetGTEX_phenotype.txt.gz",
        expand(ROOT + "data/regulators/{f}", f=["Homo_sapiens_TF.txt", "Homo_sapiens_Cof.txt", "DatabaseExtract_v_1.01.txt"])
    output: VST, REGULATORS, ONLY_TFS
    shell: "Rscript {REPO}/5_MRA/5_00_expression_regulators.R"

# Step 9 (ARACNe-AP network: hours)
rule aracne:
    input: VST, REGULATORS
    output: NETWORK
    shell: "cd {MRA} && bash {REPO}/5_MRA/aracneap.sh"

# Step 10 (also builds the regulon object from the network)
rule mra:
    input: GEO_1ROW, DE_AE, DE_GTEX, VST, NETWORK, ONLY_TFS
    output:
        CORE_MRA, SIGNATURE, MRA + "mra_results.rds", MRA + "top_15n15_mrs.tsv",
        MRA + "Figure_S1.pdf", MRA + "Figure_S2.pdf", MRA + "S1_S2_concordance.tsv"
    shell: "Rscript {REPO}/5_MRA/5_0_ovary_mra.R"

# Step 11
rule elbow:
    input: ALLGENES, CORE_MRA, REGULON
    output:
        ELBOW_TAB, ELBOW_REG, MRA + "5_1_regulon_genes_pool.rds", MRA + "5_2_core_genes_filtered.rds",
        expand(MRA + "5_0_Elbow_plot_regulons.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/5_MRA/5_1_Elbow_TF.R"

# Step 12
rule depmap_consensus:
    input: ALLGENES, CRISPR
    output:
        DEP_CORE,
        expand(DEP + "{f}.{e}", f=["6_1_1_Heatmap_coregenes", "6_1_2_RankedProfile_allgenes"], e=["pdf", "png"])
    shell: "Rscript {REPO}/6_Depmap_ovary/6_1_Heatmap_coregenes.R"

rule enriched_regulators:
    input: CORE_MRA, REGULON, ELBOW_REG, DEP_CORE
    output:
        DEP + "6_2_0_names_TF_DEG_TopNES.rds", DEP + "6_2_1_regulons_sig_filtered.rds",
        expand(DEP + "6_2_1_Enriched_TF_Depndency.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/6_Depmap_ovary/6_2_MR_Enriched.R"

rule candidate_regulators:
    input:
        DEP + "6_2_1_regulons_sig_filtered.rds", DEP + "6_2_0_names_TF_DEG_TopNES.rds", CRISPR, DEPMAP_TPM
    output: CANDIDATES, expand(DEP + "6_3_1_Heatmap_TF.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/6_Depmap_ovary/6_3_Heatmap_TF.R"

# Step 13 (both scripts write to the working directory)
rule circos_regulon_overlap:
    input: NETWORK, CORE_MRA, CANDIDATES
    output: expand(MRA + "5_2_1_circos_drivers_regulon_overlap.{e}", e=["pdf", "png"])
    shell: "cd {MRA} && Rscript {REPO}/5_MRA/5_2_Circos_regulon_overlap.R"

rule circos_driver_hallmark:
    input: NETWORK, CANDIDATES
    output: expand(MRA + "5_3_1_circos_drivers_hallmark.{e}", e=["pdf", "png"])
    shell: "cd {MRA} && Rscript {REPO}/5_MRA/5_3_Circos_driver_hallmark.R"

# Step 14 (the caches that 7_1 and 7_3 keep for themselves are not declared; 7_3 and 8_1 read those of 7_2)
rule figure2:
    input: DE_GTEX, DE_AE, GEO_RRA, ALLGENES, CONCORDANT, FANO, CNT + "3_3_2_dds_tumor_Ctlhomolog.rds", ORA
    output:
        expand(FIG + "7_1_Figure2_consensus_signature.{e}", e=["pdf", "png"]), FIG + "7_1_Figure2_panelD_heatmap.png"
    shell: "Rscript {REPO}/7_Figures/7_1_Figure2_consensus_signature.R"

rule figure3:
    input: ALLGENES, CORE_MRA, REGULON, ELBOW_REG, CRISPR, DEPMAP_TPM
    output:
        MRA_CACHE,
        expand(FIG + "7_2_Figure3_depmap_dependency.{e}", e=["pdf", "png"]), FIG + "7_2_Figure3_panelC_heatmap.png"
    shell: "Rscript {REPO}/7_Figures/7_2_Figure3_depmap_dependency.R"

rule figure4:
    input: REGULON, DEP_CORE, DEP_CACHE, MRA_CACHE
    output: expand(FIG + "7_3_Figure4_circos_master_regulators.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/7_Figures/7_3_Figure4_circos_master_regulators.R"

# Step 15 (8_0 downloads the clinical supplements from the GDC)
rule clinical_metadata:
    output: SURV + "8_0_clinical_TCGA_OV.rds", SURV + "8_0_clinical_TCGA_OV.tsv"
    shell: "Rscript {REPO}/8_Survival_analysis/8_0_Clinical_metadata.R"

rule cox_tcga:
    input:
        SURV + "8_0_clinical_TCGA_OV.rds", CNT + "3_3_2_norm_counts_tumor_ONLY.rds", CNT + "3_3_1_metadata_all.rds",
        CNT + "3_3_2_gene_annotation_UNION_universe.rds", DEP_CACHE, MRA_CACHE, CORE_MRA
    output: COX_TCGA, expand(SURV + "8_1_1_forest_cox_drivers.{e}", e=["pdf", "png"])
    shell: "Rscript {REPO}/8_Survival_analysis/8_1_Cox_drivers_survival.R"


# ---------------------------------------------------------------------------
# Revision analyses (9_Revision_analyses/)
# ---------------------------------------------------------------------------

rule rev_fano_tests:
    input: FANO, ALLGENES
    output: REV + "consensus/fano_tests.tsv", REV + "consensus/fano_paired_cohorts.tsv"
    shell: "Rscript {RSRC}/9_1_1_Fano_tests.R"

rule rev_ora_by_direction:
    input: ALLGENES
    output: REV + "consensus/ORA_GO_consensus_by_direction.tsv"
    shell: "Rscript {RSRC}/9_1_2_ORA_by_direction.R"

rule rev_prioritization_sensitivity:
    input: CORE_MRA, REGULON, ELBOW_TAB, DEP_CORE, DEP_CACHE, CANDIDATES
    output:
        expand(REV + "prioritization/{f}.tsv",
               f=["fisher_by_regulator", "HD_targets_by_regulator", "prioritization_sensitivity",
                  "null_random_HD_consensus", "null_random_HD_consensus_up"])
    shell: "Rscript {RSRC}/9_2_1_Prioritization_sensitivity.R"

rule rev_lineage_selectivity:
    input: CRISPR, CANDIDATES
    output:
        REV + "prioritization/lineage_selectivity_candidates.tsv", REV + "prioritization/ovarian_cell_lines.tsv",
        expand(REV + "prioritization/Figure_S_lineage.{e}", e=["pdf", "png"])
    shell: "Rscript {RSRC}/9_2_2_Lineage_selectivity.R"

rule rev_export_regulons:
    input: REGULON, CANDIDATES
    output: REV + "perturbseq/enriched_regulons.tsv"
    shell: "Rscript {RSRC}/9_3_1_Export_regulons.R"

rule rev_perturbseq:
    input:
        REV + "perturbseq/enriched_regulons.tsv",
        expand(REPLOGLE + "{d}_normalized_bulk_01.h5ad", d=["K562_essential", "K562_gwps", "rpe1"])
    output: REV + "perturbseq/perturbseq_regulon_tests.tsv"
    shell: "python {RSRC}/9_3_2_Perturbseq_regulon_response.py"

rule rev_perturbseq_figure:
    input: REV + "perturbseq/perturbseq_regulon_tests.tsv"
    output: expand(REV + "perturbseq/Figure_S_perturbseq.{e}", e=["pdf", "png"])
    shell: "Rscript {RSRC}/9_3_3_Figure_perturbseq.R"

# Downloads the candidate GEO series
rule rev_ft_series_metadata:
    output:
        expand(REV + "fallopian_tube/geo_data/{g}_series_matrix.txt.gz", g=["GSE10971", "GSE69428"]),
        REV + "fallopian_tube/pheno_candidate_series.rds"
    shell: "Rscript {RSRC}/9_4_1_FT_series_metadata.R"

rule rev_ft_limma:
    input: expand(REV + "fallopian_tube/geo_data/{g}_series_matrix.txt.gz", g=["GSE10971", "GSE69428"])
    output:
        REV + "fallopian_tube/de_all.rds", REV + "fallopian_tube/de_Stouffer.tsv",
        REV + "fallopian_tube/samples_by_series.tsv", REV + "fallopian_tube/series_used.tsv"
    shell: "Rscript {RSRC}/9_4_2_FT_limma.R"

rule rev_ft_consensus_concordance:
    input: ALLGENES, REV + "fallopian_tube/de_all.rds", REV + "fallopian_tube/de_Stouffer.tsv"
    output:
        REV + "fallopian_tube/consensus_concordance.tsv", REV + "fallopian_tube/consensus_concordance_by_gene.tsv"
    shell: "Rscript {RSRC}/9_4_3_FT_consensus_concordance.R"

rule rev_ft_mra:
    input:
        REGULON, CORE_MRA, SIGNATURE, ELBOW_TAB, CANDIDATES,
        REV + "fallopian_tube/de_all.rds", REV + "fallopian_tube/de_Stouffer.tsv"
    output:
        expand(REV + "fallopian_tube/{f}", f=["ft_mra_results.rds", "ft_mra_vs_main.tsv", "ft_mra_candidates.tsv",
                                               "ft_mra_summary_by_contrast.tsv"])
    shell: "Rscript {RSRC}/9_4_4_FT_mra.R"

rule rev_ft_activated_vs_repressed:
    input:
        REV + "fallopian_tube/ft_mra_vs_main.tsv", REV + "fallopian_tube/consensus_concordance.tsv",
        REV + "fallopian_tube/de_all.rds", REV + "fallopian_tube/de_Stouffer.tsv", SIGNATURE
    output:
        REV + "fallopian_tube/elbow_activated_vs_repressed.tsv", REV + "fallopian_tube/ft_mra_vs_main_ranked.tsv"
    shell: "Rscript {RSRC}/9_4_5_FT_activated_vs_repressed.R"

rule rev_ft_figure:
    input: REV + "fallopian_tube/ft_mra_vs_main.tsv", REV + "fallopian_tube/ft_mra_summary_by_contrast.tsv"
    output: expand(REV + "fallopian_tube/Figure_S_NES_main_vs_fallopian_tube.{e}", e=["pdf", "png"])
    shell: "Rscript {RSRC}/9_4_6_FT_figure_NES.R"

rule rev_cox_meta_independent:
    input: CANDIDATES, COX_TCGA
    output:
        REV + "survival_independent/results.rds",
        expand(REV + "survival_independent/{f}.tsv",
               f=["datasets_used", "cox_by_dataset", "meta_by_gene", "meta_sensitivity", "leave_one_out",
                  "group_tests", "comparison_TCGA", "aliases_used"])
    shell: "Rscript {RSRC}/9_5_1_Cox_meta_independent.R"

rule rev_forest_independent:
    input: REV + "survival_independent/results.rds"
    output:
        expand(REV + "survival_independent/{f}.{e}",
               f=["forest_meta_independent", "forest_meta_independent_high_grade"], e=["pdf", "png"])
    shell: "Rscript {RSRC}/9_5_2_Forest_independent.R"

# Needs the octad package, as steps 1, 2 and 4
rule rev_octad_reference_metric:
    output:
        expand(REV + "octad_metric/{f}.tsv",
               f=["metric_by_sample", "metric_by_tissue", "group_summary", "composition_AE_controls",
                  "composition_AE_controls_by_group", "tumor_vs_tumor_context"])
    shell: "Rscript {RSRC}/9_6_1_Octad_reference_metric.R"

rule rev_octad_reference_figure:
    input: REV + "octad_metric/metric_by_sample.tsv", REV + "octad_metric/tumor_vs_tumor_context.tsv"
    output: expand(REV + "octad_metric/octad_reference_metric.{e}", e=["pdf", "png"])
    shell: "Rscript {RSRC}/9_6_2_Octad_reference_figure.R"

# One XML per gene of 9_7_0_HPA_genes.tsv; files already present are not downloaded again
rule rev_hpa_download:
    output: expand(REV + "hpa/xml/{g}.xml", g=HPA_GENES)
    shell: "bash {RSRC}/9_7_1_HPA_download_xml.sh"

rule rev_hpa_extract:
    input: expand(REV + "hpa/xml/{g}.xml", g=HPA_GENES)
    output: REV + "hpa/hpa_v25_ovary_by_antibody.tsv", REV + "hpa/hpa_v25_ovary_patients.tsv"
    shell: "python {RSRC}/9_7_2_HPA_extract.py"

rule rev_hpa_pattern_rule:
    input: REV + "hpa/hpa_v25_ovary_by_antibody.tsv", REV + "hpa/hpa_v25_ovary_patients.tsv", CANDIDATES
    output: REV + "hpa/pattern_rule_by_antibody.tsv"
    shell: "python {RSRC}/9_7_3_HPA_pattern_rule.py"


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

PIPELINE_RULES = [r for r in workflow.rules if not r.name.startswith("rev_")]
REVISION_RULES = [r for r in workflow.rules if r.name.startswith("rev_")]

rule pipeline:
    input: [r.output for r in PIPELINE_RULES]

rule revision:
    input: [r.output for r in REVISION_RULES]

rule all:
    default_target: True
    input: rules.pipeline.input, rules.revision.input
