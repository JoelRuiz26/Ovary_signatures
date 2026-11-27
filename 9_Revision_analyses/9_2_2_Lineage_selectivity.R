#!/usr/bin/env Rscript

# ============================================================
# Lineage selectivity of the candidate regulators in DepMap, and table
# of the ovarian cell lines.
#
# For each candidate, compares the CRISPR gene effect of the ovarian
# epithelial tumor lines (and of the high-grade serous subset) with all
# other lines, and classifies the dependency by its frequency across lines.
#
# Input:  DepMap Model.csv and CRISPRGeneEffect.csv (DEPMAP_DIR),
#         6_Depmap_ovary/6_3_2_candidate_regulators.tsv
# Output (9_Revision_analyses/prioritization/):
#         lineage_selectivity_candidates.tsv, ovarian_cell_lines.tsv,
#         Figure_S_lineage.pdf / .png
# ============================================================

suppressMessages({
  library(data.table)
  library(ggplot2)
})
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

BASE       <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
DEPMAP_DIR <- Sys.getenv("DEPMAP_DIR", file.path(BASE, "data/depmap_24Q4"))
OUT_DIR    <- file.path(BASE, "9_Revision_analyses/prioritization")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

model <- fread(file.path(DEPMAP_DIR, "Model.csv"))
ge    <- fread(file.path(DEPMAP_DIR, "CRISPRGeneEffect.csv"))
setnames(ge, 1, "ModelID")
gm <- as.matrix(ge[, -1])
rownames(gm) <- ge$ModelID
colnames(gm) <- toupper(sub(" \\(\\d+\\)$", "", colnames(gm)))

ovm <- model[OncotreePrimaryDisease == "Ovarian Epithelial Tumor" & ModelID %in% rownames(gm)]
ov  <- rownames(gm) %in% ovm$ModelID
hg  <- rownames(gm) %in% ovm[OncotreeSubtype == "High-Grade Serous Ovarian Cancer"]$ModelID
cat("lines with CRISPR data:", nrow(gm), "| ovarian:", sum(ov), "| HGSOC:", sum(hg), "| other:", sum(!ov), "\n")
print(ovm[, .N, by = OncotreeSubtype][order(-N)])
cat("lineages among the other lines:", uniqueN(model[ModelID %in% rownames(gm)[!ov]]$OncotreeLineage), "\n")

cand <- fread(file.path(BASE, "6_Depmap_ovary/6_3_2_candidate_regulators.tsv"))[candidate == TRUE]$gene

tab <- rbindlist(lapply(cand, function(g) {
  x <- gm[ov, g]; y <- gm[!ov, g]; h <- gm[hg, g]
  data.table(gene = g, median_ovarian = median(x, na.rm = TRUE), median_other = median(y, na.rm = TRUE),
             median_HGSOC = median(h, na.rm = TRUE),
             pct_dependent_ovarian = 100 * mean(x < -0.5, na.rm = TRUE),
             pct_dependent_all = 100 * mean(gm[, g] < -0.5, na.rm = TRUE),
             p_two_sided = wilcox.test(x, y)$p.value,
             p_stronger_in_ovarian = wilcox.test(x, y, alternative = "less")$p.value,
             p_HGSOC_two_sided = wilcox.test(h, y)$p.value,
             p_HGSOC_stronger = wilcox.test(h, y, alternative = "less")$p.value)
}))
for (v in grep("^p_", names(tab), value = TRUE)) tab[, (sub("^p_", "padj_", v)) := p.adjust(get(v), "BH")]
tab[, delta := median_ovarian - median_other]
tab[, class := fifelse(pct_dependent_all >= 90, "common essential",
                       fifelse(pct_dependent_all >= 50, "frequent dependency", "selective dependency"))]
setorder(tab, -pct_dependent_all)
print(tab[, .(gene, ov = round(median_ovarian, 2), other = round(median_other, 2), hgsoc = round(median_HGSOC, 2),
              pct_all = round(pct_dependent_all), padj2 = signif(padj_two_sided, 2),
              padj1 = signif(padj_stronger_in_ovarian, 2), padjH2 = signif(padj_HGSOC_two_sided, 2),
              padjH1 = signif(padj_HGSOC_stronger, 2), class)])
print(table(tab$class))
fwrite(tab, file.path(OUT_DIR, "lineage_selectivity_candidates.tsv"), sep = "\t")

# Table of the ovarian lines with the gene effect of each candidate
lines <- cbind(ovm[, .(ModelID, cell_line = StrippedCellLineName, OncotreeSubtype, OncotreeCode,
                       PrimaryOrMetastasis, Sex, Age)],
               round(gm[ovm$ModelID, cand], 3))
fwrite(lines, file.path(OUT_DIR, "ovarian_cell_lines.tsv"), sep = "\t")

# Figure: gene effect of each candidate, ovarian lines against all others
long <- rbindlist(lapply(cand, function(g) data.table(
  gene = g, effect = gm[, g],
  group = fifelse(hg, "High-grade serous (ovarian)", fifelse(ov, "Other ovarian epithelial", "Other lineages")))))
long <- long[!is.na(effect)]
long[, gene := factor(gene, levels = tab[order(median_ovarian)]$gene)]
long[, group := factor(group, levels = c("Other lineages", "Other ovarian epithelial", "High-grade serous (ovarian)"))]

# Arimo has the metrics of Arial; the device falls back to the default sans if it is not installed
p <- ggplot(long, aes(gene, effect, fill = group)) +
  geom_hline(yintercept = c(-0.5, -1), linetype = c("dashed", "dotted"), colour = "grey40") +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, linewidth = 0.25, width = 0.75, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("grey80", "#9ecae1", "#de2d26"), name = NULL) +
  labs(x = NULL, y = "CRISPR gene effect (Chronos, DepMap 24Q4)") +
  theme_classic(base_size = 9, base_family = "Arimo") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "top")
ggsave(file.path(OUT_DIR, "Figure_S_lineage.pdf"), p, width = 7.2, height = 3.6, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure_S_lineage.png"), p, width = 7.2, height = 3.6, dpi = 300)
