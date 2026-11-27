#!/usr/bin/env Rscript

# ============================================================
# Forest plot of the meta-analytic HR (random effects) of the candidate
# regulators in independent cohorts, with the TCGA HR as reference.
#
# Input:  9_Revision_analyses/survival_independent/results.rds (9_5_1)
# Output (9_Revision_analyses/survival_independent/):
#         forest_meta_independent.pdf / .png              (main analysis)
#         forest_meta_independent_high_grade.pdf / .png   (sensitivity, grade 2-3)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))  # render en dashes and superscripts correctly

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
out_dir <- file.path(BASE, "9_Revision_analyses/survival_independent")
res     <- readRDS(file.path(out_dir, "results.rds"))

BASE_SIZE <- 13
FONT      <- "Arimo"   # metrics of Arial; the device falls back to the default sans if it is not installed
col_indep <- "#1A1A1A"
col_tcga  <- "#D55E00"

tcga <- res$comp[, .(gene, HR = HR_TCGA, lo = HR_TCGA_lo95, hi = HR_TCGA_hi95)]

forest <- function(m, file_stub, title) {
  m <- copy(m)[order(HR_random)]
  lev <- rev(m$gene)
  m[, gene := factor(gene, levels = lev)]

  pts <- rbind(
    m[, .(gene, HR = HR_random, lo = HR_random_lo95, hi = HR_random_hi95,
          cohort = "Independent cohorts (random-effects meta-analysis)")],
    tcga[, .(gene = factor(gene, levels = lev), HR, lo, hi, cohort = "TCGA-OV (discovery cohort)")]
  )
  pts[, cohort := factor(cohort, levels = unique(cohort))]
  # TCGA is drawn half a row below so that it does not hide the independent estimate
  pts[, y := as.integer(gene) + ifelse(as.integer(cohort) == 1, 0.14, -0.20)]

  y_lim <- c(0.3, length(lev) + 1.5)
  p1 <- ggplot(pts, aes(x = HR, y = y, colour = cohort, shape = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, yend = y, linewidth = cohort)) +
    geom_point(aes(size = cohort), fill = "white", stroke = 0.9) +
    scale_colour_manual(values = c(col_indep, col_tcga), name = NULL) +
    scale_shape_manual(values = c(16, 23), name = NULL) +
    scale_size_manual(values = c(2.6, 2.1), name = NULL) +
    scale_linewidth_manual(values = c(0.7, 0.45), name = NULL) +
    scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = c(0, 0), limits = y_lim) +
    scale_x_continuous(trans = "log", breaks = c(0.7, 0.8, 0.9, 1, 1.1, 1.2)) +
    labs(x = "Hazard ratio for death per SD of expression", y = NULL, title = title) +
    theme_classic(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      axis.text.y     = element_text(face = "bold", colour = "grey10"),
      axis.title.x    = element_text(face = "bold"),
      legend.position = "top", legend.direction = "vertical",
      legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
      legend.text     = element_text(size = BASE_SIZE - 2),
      plot.title      = element_text(size = BASE_SIZE, face = "bold"),
      plot.title.position = "plot",
      plot.margin     = margin(10, 6, 10, 10)
    )

  tab <- m[, .(gene,
               hr = sprintf("%.2f (%.2f–%.2f)", HR_random, HR_random_lo95, HR_random_hi95),
               k  = as.character(k),
               i2 = sprintf("%.0f", I2),
               bh = ifelse(BH_random < 0.001, "< 0.001", sprintf("%.3f", BH_random)))]
  tab[, y := as.integer(gene)]
  head_y <- length(lev) + 1.1
  col_x <- c(hr = 0, k = 1.10, i2 = 1.36, bh = 1.80)
  p2 <- ggplot(tab, aes(y = y)) +
    geom_text(aes(x = col_x["hr"], label = hr), hjust = 0, size = 3.5, family = FONT, colour = "grey15") +
    geom_text(aes(x = col_x["k"],  label = k),  hjust = 0, size = 3.5, family = FONT, colour = "grey15") +
    geom_text(aes(x = col_x["i2"], label = i2), hjust = 0, size = 3.5, family = FONT, colour = "grey15") +
    geom_text(aes(x = col_x["bh"], label = bh), hjust = 0, size = 3.5, family = FONT, colour = "grey15") +
    annotate("text", x = col_x, y = head_y, hjust = 0, size = 3.6, fontface = "bold", family = FONT,
             label = c("HR (95% CI)", "k", "I² (%)", "BH p")) +
    scale_x_continuous(limits = c(0, 2.3)) +
    scale_y_continuous(expand = c(0, 0), limits = y_lim) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = BASE_SIZE, base_family = FONT) +
    theme(plot.margin = margin(10, 10, 10, 0))

  fig <- p1 + p2 + plot_layout(widths = c(2.3, 1.25))
  ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), fig, width = 10, height = 7.4, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(file_stub, ".png")), fig, width = 10, height = 7.4, dpi = 400,
         bg = "white", device = ragg::agg_png)
  cat("[OK]", file_stub, "\n")
}

caption <- function(m, what) sprintf("%s: %d cohorts, up to %s patients and %d deaths per gene",
                                     what, max(m$k), format(max(m$n_patients), big.mark = ","),
                                     max(m$n_deaths))

m_main <- res$sens[analysis == "main"]
m_hg   <- res$sens[analysis == "high_grade_2_3"]

forest(m_main, "forest_meta_independent", caption(m_main, "Serous ovarian carcinoma, overall survival"))
forest(m_hg, "forest_meta_independent_high_grade",
       caption(m_hg, "Sensitivity analysis, grade 2–3 serous carcinoma"))
