#!/usr/bin/env Rscript

# ============================================================
# Supplementary figure: response of the regulons to the knockdown of
# their regulator (Replogle et al. 2022), raw and after removing the
# response shared by all perturbations.
#
# Input:  9_Revision_analyses/perturbseq/perturbseq_regulon_tests.tsv (9_3_2)
# Output: 9_Revision_analyses/perturbseq/Figure_S_perturbseq.pdf / .png
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/perturbseq")

d <- read_tsv(file.path(OUT_DIR, "perturbseq_regulon_tests.tsv"), show_col_types = FALSE) %>% filter(in_screen)
lab <- c(rpe1 = "RPE1 (essential-gene screen)", K562_essential = "K562 (essential-gene screen)",
         K562_gwps = "K562 (genome-wide screen)")
long <- bind_rows(
  d %>% transmute(dataset, driver, analysis = "Raw z-score", shift = med_z_pos - med_z_bg, padj = p_down_BH),
  d %>% transmute(dataset, driver, analysis = "After removing the shared response",
                  shift = med_zres_pos - med_zres_bg, padj = p_down_resid_BH)) %>%
  mutate(dataset = factor(lab[dataset], levels = lab),
         analysis = factor(analysis, levels = c("Raw z-score", "After removing the shared response")),
         sig = padj < 0.05)
ord <- long %>% filter(analysis == "Raw z-score", dataset == lab["rpe1"]) %>% arrange(shift) %>% pull(driver)
long$driver <- factor(long$driver, levels = c(ord, setdiff(unique(long$driver), ord)))
print(long %>% group_by(dataset, analysis) %>% summarise(n = n(), n_sig = sum(sig), .groups = "drop"))

# Arimo has the metrics of Arial; the device falls back to the default sans if it is not installed
p <- ggplot(long, aes(driver, shift, fill = sig)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_col(width = 0.7) +
  facet_grid(analysis ~ dataset, scales = "free_y") +
  scale_fill_manual(values = c(`FALSE` = "grey75", `TRUE` = "#de2d26"),
                    labels = c("n.s.", "BH-adjusted p < 0.05"), name = NULL) +
  labs(x = NULL, y = "Median expression change of positive regulon targets\nminus that of non-target genes (z)") +
  theme_bw(base_size = 8, base_family = "Arimo") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "italic", size = 6), legend.position = "top",
        panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
ggsave(file.path(OUT_DIR, "Figure_S_perturbseq.pdf"), p, width = 7.2, height = 4.4, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure_S_perturbseq.png"), p, width = 7.2, height = 4.4, dpi = 300)

print(d %>% filter(dataset == "rpe1") %>%
        select(driver, n_cells, kd_fold_expr, med_z_pos, med_z_bg, p_down_BH, med_zres_pos, med_zres_bg,
               p_down_resid_BH, own_rank) %>%
        as.data.frame(), digits = 3)
