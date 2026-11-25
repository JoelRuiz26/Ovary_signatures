#!/usr/bin/env Rscript

# ============================================================
# Figure 4 | Regulatory circuit of the essential master regulators.
#
# Takes the master regulators that are (i) enriched in the highly
# dependent fraction of the consensus signature and (ii) themselves
# essential in ovarian cancer cell lines (mean CRISPR dependency
# <= -0.5), and draws a chord diagram linking each of them to the
# highly dependent consensus genes it regulates.
#
# Inputs are the caches written by 7_2_Figure3_depmap_dependency.R
# plus the ARACNe regulons; the extracted network is cached here as
# well, so redrawing the figure takes seconds.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(vroom)
  library(purrr)
  library(circlize)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
out_dir  <- file.path(base_dir, "7_Figures")

paths <- list(
  regulons = file.path(base_dir, "5_MRA/cancer_ovary_regulon_300bt_p1e-8.rds"),
  dep_deg  = file.path(base_dir, "6_Depmap_ovary/6_1_0_Dependency_DEG_core.tsv")
)

dep_cache <- file.path(out_dir, "7_2_depmap_ovary_24Q4_646_cache.rds")
mra_cache <- file.path(out_dir, "7_2_mra_enrichment_24Q4_646_cache.rds")
net_cache <- file.path(out_dir, "7_3_circos_network_cache.rds")

DEP_MR     <- -0.5   # essentiality threshold applied to the regulators
DEP_TARGET <- -1     # highly dependent consensus genes used as targets

FIG_W <- 14.2        # inches
FIG_H <- 11.0

# The canvas is widened to the right of the circle (which always
# lives inside a radius of 1) to leave room for the legends
CANVAS_X <- c(-1, 1.58)
CANVAS_Y <- c(-1, 1)

# Dependency colour code, identical to Figure 3
col_dep_deep <- "#67001F"
col_dep_neg  <- "#B2182B"
col_dep_mid  <- "#F7F7F7"
col_dep_pos  <- "#2166AC"

col_target   <- "#C9CCD1"   # sector colour of the target genes
col_lfc      <- "#3E5C76"   # bars of the expression fold change

# ============================================================
# 1) Network: essential master regulators -> dependent targets
# ============================================================
if (file.exists(net_cache) && file.mtime(net_cache) > file.mtime(mra_cache)) {
  net <- readRDS(net_cache)
} else {
  dep <- readRDS(dep_cache)
  mra <- readRDS(mra_cache)

  # Mean dependency of every enriched regulator across the cell lines
  mr_dep <- dep$crispr %>%
    filter(gene %in% mra$results$regulon) %>%
    group_by(gene) %>%
    summarise(mean_dep = mean(dep, na.rm = TRUE), .groups = "drop") %>%
    filter(mean_dep <= DEP_MR) %>%
    arrange(mean_dep)

  # Consensus genes that are highly dependent, with their fold change
  targets <- vroom(path.expand(paths$dep_deg), show_col_types = FALSE) %>%
    filter(median_dep <= DEP_TARGET) %>%
    select(gene, median_dep, median_log2FC) %>%
    distinct()

  regulons <- readRDS(path.expand(paths$regulons))
  regulons <- regulons[mr_dep$gene]

  links <- imap_dfr(regulons, function(reg, tf) {
    tibble(
      from = tf,
      to   = intersect(names(reg$tfmode), targets$gene)
    )
  }) %>%
    filter(from != to) %>%
    mutate(value = 1)

  net <- list(mr = mr_dep, targets = targets, links = links)

  saveRDS(net, net_cache)
  rm(regulons); gc()
}

mr_genes     <- net$mr$gene
target_genes <- setdiff(unique(net$links$to), mr_genes)

cat("[INFO] essential master regulators:", length(mr_genes),
    "| dependent targets:", length(target_genes),
    "| regulatory links:", nrow(net$links), "\n")

# ============================================================
# 2) Sector order, colours and per-gene annotation values
# ============================================================
# Regulators first, from the most essential; targets afterwards,
# also sorted by dependency, so the circle reads from strong to
# weak essentiality on both arcs.
target_order <- net$targets %>%
  filter(gene %in% target_genes) %>%
  arrange(median_dep) %>%
  pull(gene)

sector_order <- c(mr_genes, target_order)

# Dependency shown by the ring: the regulators contribute their mean
# across cell lines, the targets their median in the consensus table
dep_value <- c(
  setNames(net$mr$mean_dep, net$mr$gene),
  setNames(net$targets$median_dep, net$targets$gene)[target_order]
)

lfc_value <- setNames(net$targets$median_log2FC, net$targets$gene)
lfc_value <- lfc_value[sector_order]
names(lfc_value) <- sector_order
lfc_max <- max(lfc_value, na.rm = TRUE)

col_dep <- colorRamp2(
  c(min(dep_value), min(dep_value) / 2, 0, 0.2),
  c(col_dep_deep, col_dep_neg, col_dep_mid, col_dep_pos)
)

# One hue per regulator: the links inherit it, so every ribbon can be
# traced back to the regulator it comes from
mr_col <- setNames(hcl.colors(length(mr_genes), palette = "Dark 3"), mr_genes)

grid_col <- c(mr_col, setNames(rep(col_target, length(target_genes)), target_genes))

link_col <- add_transparency(mr_col[net$links$from], 0.55)

# Gaps: wide between the two arcs, and wide enough between the narrow
# target sectors to keep their labels apart. Both groups carry the
# same number of link ends, so the regulator gaps are scaled to make
# the two arcs span the same angle and the split fall on the vertical
# axis of the figure.
gap_target <- 1.15
gap_mr     <- gap_target * (length(target_genes) - 1) / (length(mr_genes) - 1)

gap_after <- c(
  rep(gap_mr, length(mr_genes) - 1), 7,
  rep(gap_target, length(target_genes) - 1), 7
)

# ============================================================
# 3) Chord diagram
# ============================================================
draw_circos <- function() {
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.after    = gap_after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.004, 0.004),
    canvas.xlim  = CANVAS_X,
    canvas.ylim  = CANVAS_Y,
    points.overflow.warning = FALSE
  )

  chordDiagram(
    net$links,
    order          = sector_order,
    grid.col       = grid_col,
    col            = link_col,
    transparency   = 0,
    directional    = 1,
    direction.type = "diffHeight",
    diffHeight     = mm_h(1.4),
    link.sort      = TRUE,
    link.decreasing = TRUE,
    link.largest.ontop = FALSE,
    annotationTrack = "grid",
    annotationTrackHeight = mm_h(3.5),
    preAllocateTracks = list(
      list(track.height = mm_h(17)),   # gene names
      list(track.height = mm_h(3.5)),  # dependency ring
      list(track.height = mm_h(7))     # fold-change bars
    )
  )

  # --- gene names -------------------------------------------------
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    s <- CELL_META$sector.index
    circos.text(
      CELL_META$xcenter, CELL_META$ylim[1], s,
      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
      cex = if (s %in% mr_genes) 0.80 else 0.62,
      font = if (s %in% mr_genes) 2 else 1,
      col  = if (s %in% mr_genes) "grey5" else "grey25"
    )
  })

  # --- dependency ring --------------------------------------------
  circos.track(track.index = 2, bg.border = NA, panel.fun = function(x, y) {
    s <- CELL_META$sector.index
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                col = col_dep(dep_value[s]), border = NA)
  })

  # --- expression fold change --------------------------------------
  circos.track(track.index = 3, bg.border = NA, panel.fun = function(x, y) {
    s <- CELL_META$sector.index
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1,
                col = "grey96", border = NA)
    v <- lfc_value[s]
    if (!is.na(v)) {
      w <- CELL_META$xrange * 0.12
      circos.rect(CELL_META$xlim[1] + w, 0, CELL_META$xlim[2] - w, v / lfc_max,
                  col = col_lfc, border = NA)
    }
  })
}

# ============================================================
# 4) Legends
# ============================================================
# Drawn in the coordinates of the circos canvas, to the right of the
# circle, so everything stays on a single device
LG_X   <- 1.07                 # left edge of the colour patches
LG_W   <- 0.075                # their width
LG_TXT <- LG_X + LG_W + 0.035  # left edge of the labels

lg_title <- function(y, label) {
  text(LG_X, y, label, adj = c(0, 0.5), cex = 1.15, font = 2, xpd = NA)
}

lg_patch <- function(y, fill, label, h = 0.042) {
  rect(LG_X, y - h, LG_X + LG_W, y + h, col = fill, border = "grey35", lwd = 0.7,
       xpd = NA)
  text(LG_TXT, y, label, adj = c(0, 0.5), cex = 1.02, xpd = NA)
}

draw_legends <- function() {
  # --- what a sector is -------------------------------------------
  lg_title(0.86, "Sector")
  lg_patch(0.775, "grey35",   "Master regulator")
  lg_patch(0.685, col_target, "Consensus target")

  # --- dependency colour bar ---------------------------------------
  lg_title(0.55, "Dependency score")
  y_bot <- 0.06
  y_top <- 0.47
  v_lo  <- min(dep_value)
  v_hi  <- 0.2
  v_seq <- seq(v_lo, v_hi, length.out = 200)
  y_seq <- seq(y_bot, y_top, length.out = 201)

  rect(LG_X, y_seq[-201], LG_X + LG_W, y_seq[-1],
       col = col_dep(v_seq), border = NA, xpd = NA)
  rect(LG_X, y_bot, LG_X + LG_W, y_top, border = "grey35", lwd = 0.7, xpd = NA)

  ticks <- c(0, -0.5, -1, -1.5, -2)
  ticks <- ticks[ticks >= v_lo & ticks <= v_hi]
  y_tick <- y_bot + (ticks - v_lo) / (v_hi - v_lo) * (y_top - y_bot)
  segments(LG_X + LG_W, y_tick, LG_X + LG_W + 0.018, y_tick,
           col = "grey35", lwd = 0.7, xpd = NA)
  text(LG_TXT, y_tick, sprintf("%.1f", ticks), adj = c(0, 0.5), cex = 1.02,
       xpd = NA)

  # --- fold-change bars --------------------------------------------
  lg_title(-0.06, expression(bold("Median log"[2]*"FC")))
  lg_patch(-0.145, col_lfc, sprintf("bar height, 0 to %.1f", lfc_max))

  # --- how to read the ribbons -------------------------------------
  lg_title(-0.31, "Ribbons")
  text(LG_X, -0.40, "coloured by the regulator\nthey originate from",
       adj = c(0, 1), cex = 1.02, xpd = NA)
}

# ============================================================
# 5) Render
# ============================================================
render_figure <- function() {
  draw_circos()
  draw_legends()
  circos.clear()
}

out_pdf <- file.path(out_dir, "7_3_Figure4_circos_master_regulators.pdf")
out_png <- file.path(out_dir, "7_3_Figure4_circos_master_regulators.png")

cairo_pdf(out_pdf, width = FIG_W, height = FIG_H)
render_figure()
invisible(dev.off())

png(out_png, width = FIG_W, height = FIG_H, units = "in", res = 400, bg = "white")
render_figure()
invisible(dev.off())

message("Figure 4 saved to: ", out_dir)
