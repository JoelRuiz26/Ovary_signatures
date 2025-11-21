library(ggplot2)
library(dplyr)

# ---- Paths ----
meta_file <- "~/Ovary_signatures/1_DGE_signature_auto/1_1_Output_rds/Meta_DE_OVARY_ALL_auto.rds"
out_png   <- "~/Ovary_signatures/1_DGE_signature_auto/1_1_Output_rds/Volcano_Meta_OVARY_auto.png"

# Umbrales
lfc_thr  <- 1
padj_thr <- 0.01

# ---- Load meta-signature ----
df <- readRDS(meta_file) %>%
  mutate(
    mlog10p = -log10(padj),
    sig = case_when(
      log2FoldChange >=  lfc_thr & padj < padj_thr ~ "up",
      log2FoldChange <= -lfc_thr & padj < padj_thr ~ "down",
      TRUE ~ "ns"
    )
  )

# ---- Volcano plot ----
p <- ggplot(df, aes(log2FoldChange, mlog10p)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1) +
  scale_color_manual(values = c(
    up   = "#B40426",   # rojo
    down = "#3B4CC0",   # azul
    ns   = "grey70"
  )) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
  labs(
    title = "Volcano plot — Meta-firma OVARY (consensus)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw(base_size = 14)

# ---- Save PNG ----
ggsave(out_png, p, width = 7, height = 5, dpi = 600)

