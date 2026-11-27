#!/usr/bin/env Rscript

# ============================================================
# Figure 2 | Consensus gene expression signature across ovarian
#            cancer transcriptomic datasets.
#
# Rebuilds the five panels from the objects already produced by
# the analysis scripts and assembles them with patchwork:
#
#   A) Venn of significant genes        <- 3_1_DGE_list_shared.R
#   B) Fano factor of the consensus set <- 3_3_Fano_poison.R
#   C) Volcano plots (GTEx / AE)        <- 3_2_Volcano_plots_core_signature.R
#   D) Heatmap of the consensus set     <- 3_2_Heatmap.R
#   E) GO over-representation analysis  <- 4_ORA_core_bar.R
#
# Panels keep the numbers of the original figures; only the
# layout, sizes and typography are re-worked here.
# ============================================================

invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))  # render Greek letters and minus signs correctly
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(vroom)
  library(stringr)
  library(scales)
  library(ggplot2)
  library(ggvenn)
  library(ggrepel)
  library(ggh4x)
  library(pheatmap)
  library(png)
  library(grid)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root

paths <- list(
  raw_deseq2 = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae_deseq2  = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  geo_rra    = file.path(base_dir, "2_DEG_GEO",  "GEO_ovarian_cancer_RRA_results.rds"),
  all_genes  = file.path(base_dir, "3_Consensus_DGE_analysis/3_1_Core_signature/3_0_2_allgenes.tsv"),
  core       = file.path(base_dir, "3_Consensus_DGE_analysis/3_1_Core_signature/3_0_1_genes_concordant.tsv"),
  fano_gtex  = file.path(base_dir, "3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_3_fano_GTEX_raw.rds"),
  fano_ae    = file.path(base_dir, "3_Consensus_DGE_analysis/3_3_Fano_Factor/3_3_3_fano_AE_raw.rds"),
  dds_gtex   = file.path(base_dir, "3_Consensus_DGE_analysis/3_0_Counts_metadata/3_3_2_dds_tumor_Ctlhomolog.rds"),
  ora        = file.path(base_dir, "4_ORA_clusterProfiler/ORA_GO_consensus3of3_ALL_BP_MF_CC_SIMPLIFIED_0.3.tsv")
)

out_dir <- file.path(base_dir, "7_Figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# VST of the consensus genes is cached: recomputing it from the dds
# object takes several minutes and is identical on every run.
vst_cache <- file.path(out_dir, "7_1_vst_core_GTEX.rds")

# Significance thresholds (same as in the upstream scripts)
P_CUTOFF   <- 0.01
LFC_CUTOFF <- 1

# Shared typography / final canvas (landscape figure)
BASE_SIZE <- 17
TAG_SIZE  <- 26
FIG_W     <- 18.5                # inches
FIG_H     <- 11.5                # inches
ROW_H     <- c(1.05, 1.10)       # relative height of the two rows

# Panel arrangement: A, B and C on the top row, D and E on the bottom
# one. Each panel gets a cell whose proportions match its own shape
# (square Venn, wide heatmap, tall GO bars).
FIG_DESIGN <- "
AABBCCCC
DDDDDEEE
"

# Column weights of the design grid. Panels A and B carry no long axis
# labels, so they need a larger share to end up visually balanced.
COL_W <- c(2.50, 2.50, 1.55, 1.55, 1.50, 1.50, 1.50, 1.50)

# Panel D is drawn to its own raster at print resolution and then
# placed in the layout. The raster is stretched to its cell, so its
# size matches the cell proportions to avoid distorting the heatmap.
HM_W <- 9.8
HM_H <- 5.0

# Shared colour code, colour-blind safe (Okabe-Ito family).
# Data source: blue = GTEx control, amber = AE control, green = GEO.
col_gtex <- "#4878A8"
col_ae   <- "#E3A63B"
col_geo  <- "#3E9C6D"
# Lighter tints of the same hues for large filled areas (facet strips)
tint_gtex <- "#BCD2E6"
tint_ae   <- "#F6E0B0"
# Direction of regulation and diverging scale of the heatmap
col_up    <- "#B2182B"
col_down  <- "#2166AC"
col_mid   <- "#F7F7F7"
# Sample groups of the heatmap annotation
col_ctrl  <- "#79C6A5"
col_tumor <- "#D9A05B"

theme_panel <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title      = element_text(face = "bold"),
      axis.text       = element_text(colour = "grey15"),
      legend.title    = element_text(face = "bold"),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.margin     = margin(12, 16, 12, 16)
    )
}

read_rds_safe <- function(p) {
  p <- path.expand(p)
  if (!file.exists(p)) stop("File not found: ", p)
  readRDS(p)
}

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

# ============================================================
# 1) PANEL A | Venn of significant genes
# ============================================================
DE_raw <- read_rds_safe(paths$raw_deseq2)
DE_ae  <- read_rds_safe(paths$ae_deseq2)
GEO    <- read_rds_safe(paths$geo_rra)

deseq_sig <- function(df) {
  df %>%
    mutate(gene = clean_gene(Symbol)) %>%
    filter(
      !is.na(gene),
      !is.na(pvalue), pvalue < P_CUTOFF,
      !is.na(log2FoldChange), abs(log2FoldChange) > LFC_CUTOFF
    ) %>%
    pull(gene) %>%
    unique() %>%
    sort()
}

geo_sig <- GEO %>%
  mutate(gene = clean_gene(gene)) %>%
  filter(!is.na(gene), !is.na(p_raw), p_raw < P_CUTOFF) %>%
  pull(gene) %>%
  unique() %>%
  sort()

venn_list <- list(
  "GTEx\ncontrol"    = deseq_sig(DE_raw),
  "AE\ncontrol"      = deseq_sig(DE_ae),
  "GEO\n(6 datasets)" = geo_sig
)

p_venn <- ggvenn(
  venn_list,
  fill_color      = c(col_gtex, col_ae, col_geo),
  fill_alpha      = 0.55,
  stroke_size     = 0.8,
  stroke_color    = "grey25",
  set_name_size   = 6.0,
  text_color      = "transparent",   # counts are drawn below
  show_percentage = FALSE
) +
  coord_fixed(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.06)) +
  theme_void(base_size = BASE_SIZE) +
  theme(plot.margin = margin(14, 22, 14, 14))

# ggvenn writes every count with the same style and places the four-digit
# ones right on the circle strokes. They are drawn here instead.
# ggvenn lays the three sets out as unit circles centred at
# (-2/3, 0.622), (2/3, 0.622) and (0, -0.622), so every label can be
# positioned in data units, inside its own region and clear of the arcs.
venn_counts <- local({
  a <- venn_list[[1]]; b <- venn_list[[2]]; d <- venn_list[[3]]
  c(
    A   = length(setdiff(a, union(b, d))),
    B   = length(setdiff(b, union(a, d))),
    C   = length(setdiff(d, union(a, b))),
    AB  = length(setdiff(intersect(a, b), d)),
    AC  = length(setdiff(intersect(a, d), b)),
    BC  = length(setdiff(intersect(b, d), a)),
    ABC = length(Reduce(intersect, list(a, b, d)))
  )
})

venn_lab <- tibble(
  region = c("A",   "B",   "C",    "AB",  "AC",   "BC",   "ABC"),
  x      = c(-1.00,  1.00,  0.00,  0.00, -0.45,   0.45,   0.00),
  y      = c( 0.62,  0.62, -1.05,  0.75, -0.05,  -0.05,   0.15),
  face   = c(rep("plain", 6), "bold"),
  size   = c(rep(4.2, 6), 4.4)
) %>%
  mutate(n = as.integer(venn_counts[region]))

p_venn <- p_venn +
  geom_text(
    data = venn_lab,
    aes(x = x, y = y, label = n, fontface = face, size = size),
    inherit.aes = FALSE, colour = "grey10", show.legend = FALSE
  ) +
  scale_size_identity()

# The Venn has no axes to align with the other panels: wrapping it keeps
# it out of the alignment grid, so it uses the whole width of its cell.
p_venn <- wrap_elements(full = p_venn)

# ============================================================
# 2) PANEL B | Fano factor of the consensus gene set
# ============================================================
fano_GTEX <- read_rds_safe(paths$fano_gtex)
fano_AE   <- read_rds_safe(paths$fano_ae)

all_genes <- vroom(path.expand(paths$all_genes), show_col_types = FALSE) %>%
  mutate(across(c(ae, gtex, geo), as.integer))

# Subset label from the presence/absence flags of the three sources
map_src <- all_genes %>%
  transmute(
    ensembl = ensembl,
    dataset = case_when(
      ae == 1 & gtex == 0 & geo == 0 ~ "AE",
      ae == 0 & gtex == 1 & geo == 0 ~ "GTEX",
      ae == 0 & gtex == 0 & geo == 1 ~ "GEO",
      ae == 1 & gtex == 1 & geo == 0 ~ "GTEX_AE",
      ae == 1 & gtex == 0 & geo == 1 ~ "AE_GEO",
      ae == 0 & gtex == 1 & geo == 1 ~ "GTEX_GEO",
      ae == 1 & gtex == 1 & geo == 1 ~ "GTEX_AE_GEO",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(ensembl, dataset)

fano_all <- bind_rows(
  fano_GTEX %>% mutate(list = "GTEx"),
  fano_AE   %>% mutate(list = "AE")
) %>%
  left_join(map_src, by = c("gene" = "ensembl")) %>%
  filter(!is.na(dataset), is.finite(fano_factor), fano_factor >= 0) %>%
  filter(
    (list == "GTEx" & grepl("GTEX", dataset)) |
      (list == "AE"  & grepl("AE",   dataset))
  )

# Same 99th percentile trimming used in 3_3_Fano_poison.R, so the
# medians (and Δ) match the stand-alone violin plot.
fano_cut  <- as.numeric(quantile(fano_all$fano_factor, probs = 0.99, na.rm = TRUE))

fano_core <- fano_all %>%
  filter(fano_factor <= fano_cut, dataset == "GTEX_AE_GEO") %>%
  mutate(list = factor(list, levels = c("GTEx", "AE")))

med_core <- fano_core %>%
  group_by(list) %>%
  summarise(med = median(fano_factor, na.rm = TRUE), .groups = "drop")

bracket_core <- tibble(
  x      = 1.5,
  tick   = 0.06,
  y_low  = min(med_core$med),
  y_high = max(med_core$med),
  delta  = max(med_core$med) - min(med_core$med)
) %>%
  mutate(label = paste0("Δ = ", round(delta, 1)))

p_fano <- ggplot(fano_core, aes(x = list, y = fano_factor, fill = list)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.60,
              linewidth = 0.4, colour = "grey25") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.4,
               fill = "white", colour = "grey25") +
  geom_segment(data = bracket_core,
               aes(x = x, xend = x, y = y_low, yend = y_high),
               inherit.aes = FALSE, linewidth = 0.7, colour = "grey30") +
  geom_segment(data = bracket_core,
               aes(x = x - tick, xend = x + tick, y = y_low, yend = y_low),
               inherit.aes = FALSE, linewidth = 0.7, colour = "grey30") +
  geom_segment(data = bracket_core,
               aes(x = x - tick, xend = x + tick, y = y_high, yend = y_high),
               inherit.aes = FALSE, linewidth = 0.7, colour = "grey30") +
  scale_fill_manual(values = c("GTEx" = col_gtex, "AE" = col_ae)) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0.01, 0.05))) +
  coord_cartesian(ylim = c(0, 12500)) +
  # in the free corridor between the two violins, above the bracket
  annotate("text", x = 1.5, y = 12200, label = bracket_core$label,
           size = 4.1, colour = "grey20") +
  labs(x = "Control set", y = "Fano factor") +
  theme_panel() +
  theme(legend.position = "none")

# ============================================================
# 3) PANEL C | Volcano plots (TCGA-OV vs GTEx and vs AE)
# ============================================================
padj_cut     <- 0.01
lfc_cut      <- 1
TOP_N        <- 5
max_neglog10 <- 301
p_floor      <- 1e-300

prep_df <- function(df) {
  padj_col <- if ("padj_safe" %in% colnames(df)) "padj_safe" else "padj"
  df %>%
    transmute(
      Gene           = as.character(Symbol),
      log2FoldChange = as.numeric(log2FoldChange),
      padj_used      = as.numeric(.data[[padj_col]])
    ) %>%
    filter(!is.na(Gene), Gene != "",
           is.finite(log2FoldChange), is.finite(padj_used), padj_used > 0) %>%
    mutate(
      padj_used = pmax(padj_used, p_floor),
      neglog10  = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10) %>%
    mutate(
      direction = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Up",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ NA_character_
      ),
      score = abs(log2FoldChange) * neglog10
    )
}

core_tbl <- vroom(path.expand(paths$core), show_col_types = FALSE) %>%
  filter(n_sources == 3) %>%
  transmute(Gene = as.character(gene)) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct()

volcano_df <- bind_rows(
  prep_df(DE_raw) %>% mutate(dataset = "GTEx control"),
  prep_df(DE_ae)  %>% mutate(dataset = "AE control")
) %>%
  mutate(
    dataset    = factor(dataset, levels = c("GTEx control", "AE control")),
    in_core    = Gene %in% core_tbl$Gene,
    plot_group = case_when(
      in_core & direction == "Up"   ~ "Core_Up",
      in_core & direction == "Down" ~ "Core_Down",
      TRUE ~ "Other"
    )
  )

# Top up- and down-regulated consensus genes labelled per contrast
label_df <- volcano_df %>%
  filter(in_core, !is.na(direction), !grepl("[-.]", Gene)) %>%
  group_by(dataset, Gene) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  group_by(dataset, direction) %>%
  slice_max(order_by = score, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

x_lim <- c(-13, 13)
x_brk <- seq(-10, 10, by = 5)

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, neglog10)) +
  geom_point(aes(colour = plot_group, alpha = plot_group), size = 0.9) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut),     linetype = "dashed", linewidth = 0.4) +
  geom_point(
    data = label_df,
    aes(colour = plot_group), size = 2, alpha = 0.95, show.legend = FALSE
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = Gene, colour = plot_group),
    size = 4.3, fontface = "bold", segment.size = 0.3,
    bg.color = "white", bg.r = 0.14,
    box.padding = 0.50, point.padding = 0.30, min.segment.length = 0.15,
    force = 2.5, force_pull = 1.5, max.iter = 30000, max.overlaps = Inf,
    nudge_x = ifelse(label_df$direction == "Up", 0.8, -0.8),
    seed = 132, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Core_Up" = col_up, "Core_Down" = col_down, "Other" = "grey78"),
    breaks = c("Core_Up", "Core_Down"),
    labels = c("Overexpressed", "Underexpressed")
  ) +
  scale_alpha_manual(
    values = c("Core_Up" = 0.9, "Core_Down" = 0.9, "Other" = 0.35),
    guide = "none"
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4.5, alpha = 1))) +
  facet_grid2(
    . ~ dataset,
    strip = strip_themed(
      background_x = elem_list_rect(fill = c(tint_gtex, tint_ae), colour = "grey40")
    )
  ) +
  scale_x_continuous(breaks = x_brk) +
  coord_cartesian(xlim = x_lim, ylim = c(0, max_neglog10 * 1.15), clip = "off") +
  labs(x = expression(bold(log[2]~"fold change")),
       y = expression(bold(-log[10]~"(adjusted p-value)")),
       colour = "Consensus genes") +
  theme_bw(base_size = BASE_SIZE) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold"),
    strip.text.x     = element_text(face = "bold", size = BASE_SIZE - 1, colour = "grey10"),
    legend.position      = "top",
    legend.justification = "right",
    legend.direction     = "horizontal",
    legend.title         = element_text(face = "bold"),
    legend.margin        = margin(0, 0, 4, 0),
    panel.spacing.x  = unit(0.9, "lines"),
    legend.box.margin = margin(t = 6),
    plot.margin      = margin(12, 16, 12, 16)
  )

# ============================================================
# 4) PANEL D | Heatmap of the consensus gene set (VST, Z-score)
# ============================================================
vst_inputs <- path.expand(c(paths$dds_gtex, paths$all_genes))
if (file.exists(vst_cache) && all(file.mtime(vst_cache) > file.mtime(vst_inputs))) {
  hm <- readRDS(vst_cache)
} else {
  suppressPackageStartupMessages(library(DESeq2))

  dds_gtex  <- read_rds_safe(paths$dds_gtex)
  core_ens  <- all_genes %>% filter(n_sources == 3) %>% pull(ensembl) %>% unique()

  vsd  <- varianceStabilizingTransformation(dds_gtex, blind = TRUE)
  mat  <- assay(vsd)
  mat  <- mat[rownames(mat) %in% core_ens, , drop = FALSE]

  anno <- data.frame(Condition = colData(dds_gtex)$condition)
  rownames(anno) <- colnames(dds_gtex)
  anno <- anno[colnames(mat), , drop = FALSE]
  anno$Condition <- factor(ifelse(grepl("^control", anno$Condition),
                                  "control", as.character(anno$Condition)))

  hm <- list(mat = mat, anno = anno)
  saveRDS(hm, vst_cache)
  rm(dds_gtex, vsd); gc()
}

# Sample dendrogram from Pearson correlation distance (as in 3_2_Heatmap.R)
d_cols <- as.dist(1 - cor(hm$mat, method = "pearson", use = "pairwise.complete.obs"))
h_cols <- hclust(d_cols, method = "average")

# pheatmap sizes its gtable in absolute units, which would distort the
# rest of the layout: the panel is written to its own raster at print
# resolution and re-imported as a grob.
hm_png <- file.path(out_dir, "7_1_Figure2_panelD_heatmap.png")

png(hm_png, width = HM_W, height = HM_H, units = "in", res = 600)
pheatmap(
  hm$mat,
  scale             = "row",
  cluster_cols      = h_cols,
  cluster_rows      = TRUE,
  treeheight_row    = 30,
  treeheight_col    = 35,
  annotation_col    = hm$anno,
  annotation_colors = list(Condition = c(control = col_ctrl, TCGA_tumor = col_tumor)),
  annotation_names_col = FALSE,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  border_color      = NA,
  fontsize          = BASE_SIZE - 1,
  color             = colorRampPalette(c(col_down, col_mid, col_up))(100)
)
invisible(dev.off())

p_heatmap <- wrap_elements(
  full = rasterGrob(readPNG(hm_png), interpolate = TRUE)
)

# ============================================================
# 5) PANEL E | GO over-representation of the consensus set
# ============================================================
parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

ora_df <- read_tsv(path.expand(paths$ora), show_col_types = FALSE) %>%
  filter(!is.na(p.adjust), !is.na(Count)) %>%
  mutate(
    Ontology      = factor(Ontology, levels = c("BP", "MF", "CC")),
    GeneRatio_num = parse_ratio(GeneRatio),
    mlog10FDR     = -log10(p.adjust),
    Description   = str_trunc(Description, width = 58)  # one line per term; full names in the ORA table
  ) %>%
  arrange(Ontology, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

p_ora <- ggplot(ora_df, aes(x = GeneRatio_num, y = Description, fill = mlog10FDR)) +
  geom_col(width = 0.72) +
  scale_fill_viridis_c(
    option = "G", direction = -1, begin = 0.25, end = 0.88,
    limits = c(2, 10), breaks = c(2, 4, 6, 8, 10),
    oob = squish, name = expression(bold(-log[10]~"(FDR)"))
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0, 0.04))) +
  facet_grid(Ontology ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Gene ratio", y = NULL) +
  theme_panel() +
  theme(
    axis.text.y       = element_text(size = BASE_SIZE - 6, colour = "grey10", lineheight = 0.85),
    axis.text.x       = element_text(size = BASE_SIZE - 2),
    strip.background  = element_rect(fill = "grey95", colour = "grey30"),
    strip.text.y.left = element_text(face = "bold", angle = 0, size = BASE_SIZE),
    strip.placement   = "outside",
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    legend.position      = "top",
    legend.justification = "right",
    legend.direction     = "horizontal",
    legend.title         = element_text(face = "bold", vjust = 0.9),
    legend.margin        = margin(0, 0, 4, 0),
    legend.key.height    = unit(0.42, "cm"),
    legend.key.width     = unit(0.7, "cm"),
    legend.text          = element_text(size = BASE_SIZE - 4)
  )

# ============================================================
# 6) ASSEMBLY
# ============================================================
figure_2 <- p_venn + p_fano + p_volcano +
  p_heatmap + p_ora +
  plot_layout(design = FIG_DESIGN, widths = COL_W, heights = ROW_H) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = TAG_SIZE, hjust = 0, vjust = 1))

out_pdf <- file.path(out_dir, "7_1_Figure2_consensus_signature.pdf")
out_png <- file.path(out_dir, "7_1_Figure2_consensus_signature.png")

ggsave(out_pdf, figure_2, width = FIG_W, height = FIG_H, units = "in", device = cairo_pdf)
ggsave(out_png, figure_2, width = FIG_W, height = FIG_H, units = "in", dpi = 400, bg = "white")

message("Figure 2 saved to: ", out_dir)
