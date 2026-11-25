#!/usr/bin/env Rscript

# ============================================================
# Figure 3 | CRISPR dependency of the consensus signature and of
#            its master regulators in ovarian cancer cell lines.
#
# Rebuilds the three panels from the objects already produced by
# the analysis scripts and assembles them with patchwork:
#
#   A) Ranked mean dependency profile  <- 6_1_Heatmap_coregenes.R
#   B) Master regulators enriched in
#      highly dependent genes          <- 6_2_MR_Enriched.R
#   C) Dependency heatmap of the
#      enriched master regulators      <- 6_3_Heatmap_TF.R
#
# Revision (FIG, 2026-11): DepMap release 24Q4 read from the
# CSV files (CRISPRGeneEffect.csv, Model.csv) instead of the
# Bioconductor 'depmap' package, which ships release 22Q2; 58
# "Ovarian Epithelial Tumor" lines; consensus = 646 genes
# significant in all three sources (direction by agreement of
# >= 2 of 3); panel A shows the
# median dependency (the statistic used for the HD-genes).
# Same logic as 6_1, 6_2 and 6_3 in branch revision-depmap24q4.
#
# The DepMap CRISPR table (~19 M rows) and the Fisher enrichment
# are cached in this folder, so re-running the script to adjust
# the figure takes seconds instead of minutes.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(vroom)
  library(purrr)
  library(forcats)
  library(scales)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(png)
  library(grid)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir   <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")   # project data root
depmap_dir <- Sys.getenv("DEPMAP_DIR", file.path(base_dir, "data/depmap_24Q4"))
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))  # render Unicode minus signs correctly

paths <- list(
  all_genes  = file.path(base_dir, "3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv"),
  mra_tsv    = file.path(base_dir, "5_MRA/core_mra.tsv"),
  regulons   = file.path(base_dir, "5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds"),
  top_nes    = file.path(base_dir, "5_MRA/5_0_Top_elbow_regulones.rds")
)

out_dir <- Sys.getenv("FIG_OUT_DIR", file.path(base_dir, "7_Figures"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Caches are tagged with the release so the old 22Q2 caches are never reused
dep_cache <- file.path(out_dir, "7_2_depmap_ovary_24Q4_646_cache.rds")
mra_cache <- file.path(out_dir, "7_2_mra_enrichment_24Q4_646_cache.rds")

# Thresholds of the upstream scripts
DEP_ESSENTIAL <- -0.5    # dependency line drawn in panel A
DEP_STRONG    <- -1      # genes taken as highly dependent for the enrichment
P_ADJ_MR      <- 0.01
MIN_OVERLAP   <- 10
TOP_N         <- 10      # genes shown in the inset of panel A

# Shared typography / final canvas (landscape figure)
BASE_SIZE <- 16
TAG_SIZE  <- 26
FIG_W     <- 18.0                    # inches
FIG_H     <- 17.5                    # inches

# Two rows: the ranked profile and the enrichment share the top row,
# and the heatmap takes the full width below so that the 58 cell
# lines and the 46 regulators are printed at a readable size.
FIG_DESIGN <- "
AAAABBB
CCCCCCC
"
COL_W <- c(rep(1.28, 4), rep(1.00, 3))
ROW_H <- c(0.95, 1.25)

# Panel C is drawn to its own raster at print resolution: its size
# matches the cell proportions to avoid distorting the heatmap.
HM_W <- 17.6
HM_H <- 9.8

# Colour code: red = essential (negative dependency), blue =
# dispensable. Panels A and C share the very same four stops, so the
# ranked profile and the heatmap can be read against each other.
col_dep_deep <- "#67001F"
col_dep_neg  <- "#B2182B"
col_dep_mid  <- "#F7F7F7"
col_dep_pos  <- "#2166AC"

# Enrichment significance: purples, the only hue not already carrying
# a meaning in the figure (red/blue = dependency, green = signature).
pal_enrich <- c("#DADAEB", "#9E9AC8", "#6A51A3", "#3F007D")

# Regulators that also belong to the consensus signature
col_dge     <- "#3E9C6D"
col_no_dge  <- "#DEDEDE"

theme_panel <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(colour = "grey15"),
      legend.title = element_text(face = "bold"),
      plot.margin  = margin(12, 16, 12, 16)
    )
}

# ============================================================
# 1) DepMap CRISPR data for ovarian cancer cell lines (cached)
# ============================================================
dep_inputs <- c(path.expand(paths$all_genes), file.path(depmap_dir, c("CRISPRGeneEffect.csv", "Model.csv")))
if (file.exists(dep_cache) && all(file.mtime(dep_cache) > file.mtime(dep_inputs))) {
  dep <- readRDS(dep_cache)
} else {
  # Significant in all three sources; direction by agreement of >= 2 of 3 (646 genes)
  cons_genes <- vroom(path.expand(paths$all_genes), show_col_types = FALSE) %>%
    filter(n_sources == 3) %>%
    transmute(gene_u = toupper(gene)) %>%
    distinct()

  # DepMap 24Q4: Chronos gene effect for the ovarian epithelial tumor lines
  ovary_lines <- vroom(file.path(depmap_dir, "Model.csv"), show_col_types = FALSE) %>%
    filter(OncotreePrimaryDisease == "Ovarian Epithelial Tumor") %>%
    dplyr::select(depmap_id = ModelID, cell = StrippedCellLineName)

  gene_effect <- vroom(file.path(depmap_dir, "CRISPRGeneEffect.csv"), show_col_types = FALSE)
  names(gene_effect)[1] <- "depmap_id"

  crispr <- gene_effect %>%
    inner_join(ovary_lines, by = "depmap_id") %>%
    dplyr::select(-depmap_id) %>%
    pivot_longer(-cell, names_to = "gene", values_to = "dep") %>%
    mutate(gene = toupper(sub(" \\(\\d+\\)$", "", gene)))
  rm(gene_effect)

  to_matrix <- function(df) {
    df %>%
      distinct() %>%
      pivot_wider(names_from = cell, values_from = dep) %>%
      column_to_rownames("gene") %>%
      as.matrix()
  }

  dep <- list(
    core      = to_matrix(filter(crispr, gene %in% cons_genes$gene_u)),
    all_genes = sort(unique(crispr$gene)),
    range     = range(crispr$dep, na.rm = TRUE),
    crispr    = crispr
  )

  saveRDS(dep, dep_cache)
  rm(crispr); gc()
}

cat("[INFO] consensus genes with dependency data:", nrow(dep$core),
    "x", ncol(dep$core), "cell lines\n")

# ============================================================
# 2) PANEL A | Ranked mean dependency profile
# ============================================================
dep_summary <- as.data.frame(dep$core) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "dep") %>%
  group_by(gene) %>%
  summarise(
    median_dep = median(dep, na.rm = TRUE),
    sd_dep     = sd(dep, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(median_dep) %>%
  mutate(rank = row_number())

dep_lims <- range(dep_summary$median_dep, na.rm = TRUE)

# Same four stops and the same anchoring as the heatmap of panel C,
# so that a given tone means the same dependency in both panels
dep_stops <- c(dep_lims[1], dep_lims[1] / 2, 0, dep_lims[2])

scale_dep <- function(aesthetic, guide) {
  do.call(
    paste0("scale_", aesthetic, "_gradientn"),
    list(
      colours = c(col_dep_deep, col_dep_neg, col_dep_mid, col_dep_pos),
      values  = rescale(dep_stops),
      limits  = dep_lims,   # same limits in the main plot and in the inset
      name    = "Median dependency score",
      guide   = guide
    )
  )
}

zoom_df <- dep_summary %>% filter(rank <= TOP_N)

p_zoom <- ggplot(zoom_df, aes(rank, median_dep)) +
  geom_errorbar(aes(ymin = median_dep - sd_dep, ymax = median_dep + sd_dep),
                width = 0.25, linewidth = 0.45, colour = "grey55") +
  geom_point(aes(fill = median_dep), shape = 21, size = 3.2,
             colour = "grey25", stroke = 0.3) +
  geom_text_repel(
    aes(label = gene), size = 3.7, fontface = "bold", colour = "grey10",
    bg.color = "white", bg.r = 0.16, segment.colour = "grey55",
    segment.size = 0.3, box.padding = 0.62, point.padding = 0.35,
    force = 3, max.overlaps = Inf, seed = 7
  ) +
  scale_dep("fill", guide = "none") +
  scale_x_continuous(breaks = seq(2, TOP_N, by = 2)) +
  labs(x = NULL, y = NULL, title = paste("Top", TOP_N, "essential genes")) +
  theme_classic(base_size = BASE_SIZE - 3) +
  theme(
    plot.title       = element_text(size = BASE_SIZE - 2, hjust = 0.5, colour = "grey20"),
    panel.background = element_rect(fill = "white", colour = "grey60"),
    plot.background  = element_rect(fill = "white", colour = NA),
    axis.line        = element_blank(),
    plot.margin      = margin(6, 8, 6, 6)
  )

p_rank <- ggplot(dep_summary, aes(rank, median_dep)) +
  geom_ribbon(aes(ymin = median_dep - sd_dep, ymax = median_dep + sd_dep, fill = median_dep),
              alpha = 0.22, colour = NA) +
  geom_line(colour = "grey45", linewidth = 0.4) +
  geom_point(aes(fill = median_dep), shape = 21, size = 1.5,
             colour = "grey25", stroke = 0.18) +
  geom_hline(yintercept = DEP_ESSENTIAL, linetype = "dotted",
             colour = col_dep_neg, linewidth = 0.8) +
  annotate("text", x = max(dep_summary$rank), y = DEP_ESSENTIAL,
           label = "dependency = -0.5", hjust = 1, vjust = -0.7,
           size = 4.8, fontface = "italic", colour = col_dep_neg) +
  scale_dep("fill", guide = guide_colourbar(title.position = "top")) +
  scale_x_continuous(breaks = c(1, seq(100, nrow(dep_summary), by = 100)),
                     expand = expansion(mult = 0.02)) +
  labs(x = "Gene rank", y = "Median dependency score across cell lines") +
  theme_panel() +
  theme(
    legend.position      = "top",
    legend.justification = "left",
    legend.direction     = "horizontal",
    legend.key.height    = unit(0.42, "cm"),
    legend.key.width     = unit(1.6, "cm"),
    legend.text          = element_text(size = BASE_SIZE - 4),
    legend.margin        = margin(0, 0, 6, 0)
  )

p_panelA <- wrap_elements(full = p_rank +
  inset_element(p_zoom, left = 0.34, bottom = 0.06, right = 1.0, top = 0.62,
                align_to = "panel"))

# ============================================================
# 3) PANEL B | Master regulators enriched in dependent genes
# ============================================================
# The cache is reused only if it is newer than the MRA table, the regulons,
# the elbow selection and the dependency cache it was computed from.
mra_inputs <- c(path.expand(unlist(paths[c("mra_tsv", "regulons", "top_nes")])), dep_cache)
if (file.exists(mra_cache) && all(file.mtime(mra_cache) > file.mtime(mra_inputs))) {
  mra <- readRDS(mra_cache)
} else {
  sig_tf <- vroom(path.expand(paths$mra_tsv), show_col_types = FALSE) %>%
    filter(padj_emp <= 0.05) %>%
    pull(TF)

  regulons <- readRDS(path.expand(paths$regulons))
  regulons <- regulons[names(regulons) %in% sig_tf]

  top_nes  <- names(readRDS(path.expand(paths$top_nes)))

  # HD-genes: consensus genes with median dependency <= -1 in the 24Q4 ovarian lines
  core_genes <- rownames(dep$core)
  genes_hi   <- dep_summary %>% filter(median_dep <= DEP_STRONG) %>% pull(gene) %>% unique()

  universe <- unique(unlist(lapply(regulons, function(x) names(x$tfmode))))
  genes_hi <- intersect(genes_hi, universe)

  # One-sided Fisher exact test per regulon
  enrich_one <- function(reg_name, regulon) {
    genes_reg <- intersect(names(regulon$tfmode), universe)
    k <- length(intersect(genes_reg, genes_hi))
    M <- length(genes_reg)
    n <- length(genes_hi)
    N <- length(universe)
    tibble(
      regulon = reg_name,
      overlap = k,
      p_value = fisher.test(matrix(c(k, M - k, n - k, N - M - n + k), nrow = 2),
                            alternative = "greater")$p.value
    )
  }

  res <- imap_dfr(regulons, ~ enrich_one(.y, .x)) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    filter(regulon %in% top_nes) %>%
    filter(p_adj < P_ADJ_MR, overlap >= MIN_OVERLAP) %>%
    mutate(log10_p = -log10(p_adj))

  mra <- list(
    results = res,
    tf_dge  = res %>% filter(regulon %in% core_genes) %>% pull(regulon),
    n_genes = length(genes_hi),
    n_hd    = sum(dep_summary$median_dep <= DEP_STRONG, na.rm = TRUE)
  )

  saveRDS(mra, mra_cache)
  rm(regulons); gc()
}

# Row labels were sized for 46 regulators; shrink them when there are more
lab_scale <- min(1, 46 / nrow(mra$results))

cat("[INFO] HD-genes (median <= ", DEP_STRONG, "): ", mra$n_hd,
    " | in the regulon universe: ", mra$n_genes, "\n", sep = "")
cat("[INFO] enriched master regulators:", nrow(mra$results),
    "| also differentially expressed:", length(mra$tf_dge), "\n")

mra_df <- mra$results %>% mutate(regulon = fct_reorder(regulon, overlap))

# Regulators that are themselves part of the consensus signature
# are highlighted in bold
mra_face <- ifelse(levels(mra_df$regulon) %in% mra$tf_dge, "bold", "plain")

p_panelB <- ggplot(mra_df, aes(x = overlap, y = regulon, fill = log10_p)) +
  geom_col(width = 0.78, colour = "grey30", linewidth = 0.15) +
  scale_fill_gradientn(colours = pal_enrich,
                       name = expression(bold(-log[10]~"(adj p)"))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Genes shared with the highly dependent set", y = NULL) +
  theme_panel() +
  theme(
    axis.text.y          = element_text(size = (BASE_SIZE - 5) * lab_scale, face = mra_face,
                                        colour = "grey10"),
    axis.text.x          = element_text(size = BASE_SIZE - 3),
    panel.grid.major.x   = element_line(colour = "grey90", linewidth = 0.3),
    legend.position      = "top",
    legend.justification = "right",
    legend.direction     = "horizontal",
    legend.key.height    = unit(0.42, "cm"),
    legend.key.width     = unit(1.1, "cm"),
    legend.text          = element_text(size = BASE_SIZE - 4),
    legend.margin        = margin(0, 0, 6, 0)
  )

# ============================================================
# 4) PANEL C | Dependency heatmap of the enriched regulators
# ============================================================
mat_tf <- dep$crispr %>%
  filter(gene %in% mra$results$regulon) %>%
  distinct() %>%
  pivot_wider(names_from = cell, values_from = dep) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Drop regulators without a 24Q4 gene effect in most lines (HMGB3:
# its multi-aligning guides were removed by DepMap)
na_frac <- rowMeans(is.na(mat_tf))
if (any(na_frac > 0.5)) cat("[INFO] dropped from panel C (mostly NA):",
                            names(na_frac)[na_frac > 0.5], "\n")
mat_tf <- mat_tf[na_frac <= 0.5, , drop = FALSE]

mat_tf <- mat_tf[order(apply(mat_tf, 1, median, na.rm = TRUE)), , drop = FALSE]  # by median, as the HD-genes

# The full CRISPR range is dominated by a few extreme values and
# flattens the map; the ramp is anchored on the 1st/99th percentiles
# of the regulators shown, which restores the contrast between
# essential and dispensable profiles.
hm_lo <- quantile(mat_tf, 0.01, na.rm = TRUE)
hm_hi <- max(mat_tf, na.rm = TRUE)

col_fun <- colorRamp2(
  c(hm_lo, hm_lo / 2, 0, hm_hi),
  c(col_dep_deep, col_dep_neg, col_dep_mid, col_dep_pos)
)

is_dge   <- rownames(mat_tf) %in% toupper(mra$tf_dge)
row_face <- ifelse(is_dge, "bold", "plain")

# Colour strip duplicating the bold labels, so the regulators that
# are also part of the consensus signature can be spotted at a glance
row_anno <- rowAnnotation(
  Signature = ifelse(is_dge, "In consensus signature", "Regulator only"),
  col = list(Signature = c(`In consensus signature` = col_dge,
                           `Regulator only`         = col_no_dge)),
  simple_anno_size   = unit(4.5, "mm"),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title     = "",
    labels_gp = gpar(fontsize = BASE_SIZE - 3),
    grid_height = unit(6, "mm"),
    grid_width  = unit(6, "mm")
  )
)

# Regulators in the rows (ordered by mean dependency, most essential
# on top) and cell lines in the columns, clustered by their profile.
# The panel spans the full width of the figure, so both label sets
# are printed large enough to be read.
ht <- Heatmap(
  mat_tf,
  name            = "Dependency\nscore",
  col             = col_fun,
  cluster_rows    = FALSE,
  cluster_columns = TRUE,
  show_row_names  = TRUE,
  row_names_gp    = gpar(fontsize = (BASE_SIZE - 2) * lab_scale, fontface = row_face),
  row_names_side  = "left",
  left_annotation = row_anno,
  show_column_names   = TRUE,
  # Vertical labels: at 45 degrees the diagonal spacing between the
  # 58 columns caps the font size, whereas at 90 degrees the whole
  # column width is available and the names can be set much larger
  column_names_gp     = gpar(fontsize = BASE_SIZE - 2, colour = "grey10"),
  column_names_rot    = 90,
  column_names_max_height = unit(55, "mm"),
  column_dend_height  = unit(14, "mm"),
  rect_gp         = gpar(col = NA),
  border          = FALSE,
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = BASE_SIZE - 1, fontface = "bold"),
    labels_gp     = gpar(fontsize = BASE_SIZE - 2),
    legend_height = unit(4.2, "cm"),
    legend_width  = unit(0.9, "cm")
  )
)

hm_png <- file.path(out_dir, "7_2_Figure3_panelC_heatmap.png")

png(hm_png, width = HM_W, height = HM_H, units = "in", res = 600)
draw(ht, heatmap_legend_side = "right", padding = unit(c(2, 2, 2, 14), "mm"))
invisible(dev.off())

p_panelC <- wrap_elements(full = rasterGrob(readPNG(hm_png), interpolate = TRUE))

# ============================================================
# 5) ASSEMBLY
# ============================================================
figure_3 <- p_panelA + p_panelB + p_panelC +
  plot_layout(design = FIG_DESIGN, widths = COL_W, heights = ROW_H) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = TAG_SIZE, hjust = 0, vjust = 1))

out_pdf <- file.path(out_dir, "7_2_Figure3_depmap_dependency.pdf")
out_png <- file.path(out_dir, "7_2_Figure3_depmap_dependency.png")

ggsave(out_pdf, figure_3, width = FIG_W, height = FIG_H, units = "in", device = cairo_pdf)
ggsave(out_png, figure_3, width = FIG_W, height = FIG_H, units = "in", dpi = 400, bg = "white")

message("Figure 3 saved to: ", out_dir)
