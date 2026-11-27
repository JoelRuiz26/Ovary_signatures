#!/usr/bin/env Rscript

# ============================================================
# Fallopian tube contrast, step 6: scatter of the NES of the main
# analysis against the NES obtained with the tumor vs fallopian tube
# signatures, highlighting the candidate regulators.
#
# Input:  9_Revision_analyses/fallopian_tube/ft_mra_vs_main.tsv,
#         ft_mra_summary_by_contrast.tsv (9_4_4)
# Output: 9_Revision_analyses/fallopian_tube/Figure_S_NES_main_vs_fallopian_tube.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
})
Sys.setlocale("LC_CTYPE", "C.UTF-8")

BASE    <- Sys.getenv("OVARY_SIGNATURES", "~/Ovary_signatures")  # project data root
OUT_DIR <- file.path(BASE, "9_Revision_analyses/fallopian_tube")

tab <- read.delim(file.path(OUT_DIR, "ft_mra_vs_main.tsv"))
res <- read.delim(file.path(OUT_DIR, "ft_mra_summary_by_contrast.tsv"))

panels <- c(GSE10971 = "GSE10971 (11 HGSC vs 12 FTE)", GSE69428 = "GSE69428 (10 paired HGSOC vs oviduct)",
            Stouffer = "Combined (Stouffer)")
lab_cand <- sprintf("%d candidates", sum(tab$candidate))
classes  <- c("Other regulators", "Above elbow, activated", "Above elbow, repressed", lab_cand)
d <- tab %>% select(TF, NES_main, above_elbow, candidate, all_of(paste0("NES_", names(panels)))) %>%
  pivot_longer(starts_with("NES_G") | starts_with("NES_S"), names_to = "contrast", values_to = "NES_ft", names_prefix = "NES_") %>%
  filter(!is.na(NES_ft)) %>%
  mutate(class = case_when(candidate ~ lab_cand,
                           above_elbow & NES_main > 0 ~ "Above elbow, activated",
                           above_elbow & NES_main < 0 ~ "Above elbow, repressed",
                           TRUE ~ "Other regulators"),
         class = factor(class, levels = classes),
         panel = factor(panels[contrast], levels = panels))
# rho per panel (all shared regulators)
rho <- res %>% filter(set == "all_common", contrast %in% names(panels)) %>%
  transmute(panel = factor(panels[contrast], levels = panels),
            label = sprintf("Spearman rho = %.2f\nn = %d regulators", rho_spearman, n_evaluated))
# the candidates are highlighted without labels (their values are in ft_mra_candidates.tsv)

# Arimo has the metrics of Arial; the device falls back to the default sans if it is not installed
cols <- setNames(c("grey75", "#D55E00", "#0072B2", "black"), classes)
p <- ggplot(d %>% arrange(class), aes(NES_main, NES_ft)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(colour = class, shape = class, size = class), fill = "#E69F00", alpha = 0.8, stroke = 0.5) +
  geom_text(data = rho, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.3, size = 2.8,
            family = "Arimo", lineheight = 0.95) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = setNames(c(16, 16, 16, 21), classes), name = NULL) +
  scale_size_manual(values = setNames(c(0.5, 1, 1, 1.8), classes), name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "NES, main analysis (tumor vs normal ovary consensus signature)",
       y = "NES, tumor vs fallopian tube epithelium signature") +
  theme_bw(base_size = 9, base_family = "Arimo") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
        legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.margin = margin(0, 0, 0, 0))
ggsave(file.path(OUT_DIR, "Figure_S_NES_main_vs_fallopian_tube.pdf"), p, width = 9.5, height = 4.1, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure_S_NES_main_vs_fallopian_tube.png"), p, width = 9.5, height = 4.1, dpi = 300)
sessionInfo()
