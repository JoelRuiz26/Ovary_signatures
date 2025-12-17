# ===================== CONFIG ===================== #
library(dplyr)
library(ggvenn)
library(ggplot2)
library(stringr)

base_dir <- "~/Ovary_signatures/"

paths <- list(
  raw_full_deseq2 = file.path(base_dir, "0_DGE_raw",         "DE_full_OVARY_DESeq2.rds"),
  ae_full_deseq2  = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds")
)

read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe el archivo: ", p)
  readRDS(p)
}

DE_full_raw_DESeq2 <- read_rds_safe(paths$raw_full_deseq2)
DE_full_ae_DESeq2  <- read_rds_safe(paths$ae_full_deseq2)

# ===================== UMBRALES ===================== #
PADJ_CUTOFF <- 0.01
LFC_CUTOFF  <- 1

# ===================== FILTRADO ===================== #
raw_sig <- DE_full_raw_DESeq2 %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF,
         !is.na(log2FoldChange), abs(log2FoldChange) > LFC_CUTOFF) %>%
  transmute(identifier)

ae_sig <- DE_full_ae_DESeq2 %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF,
         !is.na(log2FoldChange), abs(log2FoldChange) > LFC_CUTOFF) %>%
  transmute(identifier)

venn_list <- list(
  ovary_gtex_control = unique(raw_sig$identifier),
  Reference_control  = unique(ae_sig$identifier)
)

# ===================== TITULO CORTO + WRAP ===================== #
title_txt <- "Shared significant DE genes in ovarian cancer (DESeq2)"
subtitle_txt <- paste0("padj < ", PADJ_CUTOFF, " and |log2FC| > ", LFC_CUTOFF)

# ===================== PLOT (bonito) ===================== #
venn_plot <- ggvenn(
  venn_list,
  fill_color      = c("#1F77B4", "#FF7F0E"),  # azul/naranja (alto contraste)
  fill_alpha      = 0.45,
  stroke_size     = 0.9,
  set_name_size   = 6,
  text_size       = 6,
  show_percentage = TRUE,
  digits          = 1,
  label_sep       = "\n"   # número arriba, % abajo
) +
  labs(
    title = str_wrap(title_txt, width = 55),
    subtitle = subtitle_txt
  ) +
  coord_fixed() +  # evita deformación
  theme_classic(base_size = 13) +
  theme(
    axis.line  = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title    = element_text(face = "bold", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10)),
    plot.margin   = margin(t = 14, r = 18, b = 12, l = 18)
  )

# ===================== GUARDAR PDF ===================== #
out_pdf <- file.path(base_dir, "venn_DESeq2.pdf")
ggsave(out_pdf, plot = venn_plot, device = cairo_pdf, width = 7.5, height = 6.5, units = "in")

cat("PDF guardado en:", out_pdf, "\n")
