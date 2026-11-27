#!/usr/bin/env Rscript

# ============================================================
# Figure: distribution of the computeRefTissue metric and tissue
# composition of the selected reference samples.
#
# Input:  9_Revision_analyses/octad_metric/metric_by_sample.tsv,
#         tumor_vs_tumor_context.tsv (9_6_1_Octad_reference_metric.R)
# Output: 9_Revision_analyses/octad_metric/octad_reference_metric.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
Sys.setlocale("LC_CTYPE", "C.UTF-8")

BASE <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
d    <- file.path(BASE, "9_Revision_analyses/octad_metric")
tab  <- read.delim(file.path(d, "metric_by_sample.tsv"), check.names = FALSE)
n_tumors <- nrow(read.delim(file.path(d, "tumor_vs_tumor_context.tsv")))
n_ae <- sum(tab$selected_AE); n_ov <- sum(tab$ovary_GTEx); n_both <- sum(tab$selected_AE & tab$ovary_GTEx)

# Arimo has the metrics of Arial; the device falls back to the default sans if it is not installed
blue <- "#2a78d6"; orange <- "#eb6834"; grey <- "#c9c8c2"; ink <- "#0b0b0b"; ink2 <- "#52514e"
cutoff <- min(tab$median_spearman[tab$selected_AE])
xl <- range(tab$median_spearman) + c(-0.01, 0.015)
base <- theme_classic(base_size = 7, base_family = "Arimo") +
  theme(axis.line = element_line(linewidth = 0.3, colour = ink2), axis.ticks = element_line(linewidth = 0.3, colour = ink2),
        axis.text = element_text(colour = ink), plot.title = element_text(size = 7.5, face = "bold"),
        plot.title.position = "panel", plot.tag = element_text(size = 9, face = "bold"))

# a) all the candidates
pA <- ggplot(tab, aes(median_spearman)) +
  geom_histogram(binwidth = 0.01, boundary = 0, fill = grey, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = cutoff, linetype = "22", linewidth = 0.35, colour = ink) +
  annotate("text", x = cutoff - 0.006, y = Inf, vjust = 1.3, hjust = 1, size = 2.1, family = "Arimo", colour = ink,
           label = sprintf("%dth-ranked sample\nρ = %.3f", n_ae, cutoff)) +
  scale_x_continuous(limits = xl, expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = sprintf("All candidate normal and adjacent samples (n = %s)", format(nrow(tab), big.mark = ",")), x = NULL, y = "Samples") +
  base + theme(axis.text.x = element_blank())

# b) the two sets of controls
lev <- c(sprintf("GTEx ovary controls (n = %d)", n_ov), sprintf("Selected by computeRefTissue (n = %d)", n_ae))
s <- bind_rows(tab %>% filter(selected_AE) %>% mutate(group = lev[2]), tab %>% filter(ovary_GTEx) %>% mutate(group = lev[1])) %>%
  mutate(group = factor(group, lev), both = selected_AE & ovary_GTEx)
set.seed(1)  # only for the vertical jitter of the points
s$y <- as.numeric(s$group) + runif(nrow(s), -0.27, 0.27)
med <- s %>% group_by(group) %>% summarise(m = median(median_spearman), .groups = "drop") %>% mutate(yy = as.numeric(group))
pB <- ggplot(s, aes(median_spearman, y)) +
  geom_vline(xintercept = cutoff, linetype = "22", linewidth = 0.35, colour = ink) +
  geom_point(aes(colour = group), size = 0.9, alpha = 0.75, stroke = 0) +
  geom_point(data = filter(s, both), shape = 21, size = 1.9, stroke = 0.45, colour = ink, fill = NA) +
  geom_segment(data = med, aes(x = m, xend = m, y = yy - 0.38, yend = yy + 0.38), linewidth = 0.5, colour = ink) +
  geom_text(data = med, aes(x = c(m[1], xl[2] - 0.003), y = yy + 0.47, label = sprintf("median %.3f", m)), hjust = c(0.5, 1), size = 2.1, family = "Arimo", colour = ink) +
  annotate("text", x = xl[1] + 0.005, y = 2.25, hjust = 0, size = 2.1, family = "Arimo", colour = ink2,
           label = sprintf("Circled: the %d GTEx ovary samples\npresent in both control sets", n_both)) +
  scale_colour_manual(values = setNames(c(orange, blue), lev), guide = "none") +
  scale_x_continuous(limits = xl, expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:2, labels = c(sprintf("GTEx ovary\ncontrols\n(n = %d)", n_ov), sprintf("Selected by\ncomputeRefTissue\n(n = %d)", n_ae)),
                     limits = c(0.55, 2.6), expand = c(0, 0)) +
  labs(x = sprintf("Median Spearman ρ with the %d TCGA-OV tumours (64 autoencoder features)", n_tumors), y = NULL) +
  base + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# c) tissue composition of the selected samples
abr <- c("GTEx" = "GTEx", "TCGA adjacent" = "TCGA adj.", "TARGET adjacent" = "TARGET adj.")
nom <- c(KIDNEY = "Kidney", BRAIN = "Brain", PITUITARY = "Pituitary", THYROID = "Thyroid", BREAST = "Breast", "UTERUS/ENDOMETRIUM" = "Uterus / endometrium",
         "FALLOPIAN TUBE" = "Fallopian tube", OVARY = "Ovary", STOMACH = "Stomach", ESOPHAGUS = "Esophagus")
tissue_label <- function(x) ifelse(x %in% names(nom), nom[x], tools::toTitleCase(tolower(x)))
cmp <- tab %>% filter(selected_AE) %>% dplyr::count(tissue_group, source) %>% group_by(tissue_group) %>% arrange(desc(n), .by_group = TRUE) %>%
  summarise(total = sum(n), det = paste(sprintf("%s %d", abr[source], n), collapse = ", "), .groups = "drop") %>%
  arrange(total, desc(tissue_group)) %>% mutate(tissue_lab = factor(tissue_label(tissue_group), tissue_label(tissue_group)))
stopifnot(!anyNA(cmp$tissue_lab), sum(cmp$total) == n_ae)
pC <- ggplot(cmp, aes(total, tissue_lab)) +
  geom_col(fill = blue, width = 0.62) +
  geom_text(aes(label = sprintf("%d  (%s)", total, sub("^([^,]+, [^,]+), ", "\\1,\n      ", det))), hjust = 0, nudge_x = 0.7, size = 2.1, lineheight = 0.95, family = "Arimo", colour = ink) +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 30, 10), expand = c(0, 0)) +
  labs(title = sprintf("Tissue of origin of the %d selected samples", n_ae), x = "Samples", y = NULL) +
  base + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

fig <- ((pA / pB) + plot_layout(heights = c(1.15, 1))) | pC
fig <- fig + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "a")
cairo_pdf(file.path(d, "octad_reference_metric.pdf"), width = 7.2, height = 3.6, family = "Arimo"); print(fig); invisible(dev.off())
ragg::agg_png(file.path(d, "octad_reference_metric.png"), width = 7.2, height = 3.6, units = "in", res = 400); print(fig); invisible(dev.off())
cat("Figure written\n"); print(as.data.frame(cmp))
