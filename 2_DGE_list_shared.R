# ===================== CONFIG ===================== #
library(dplyr)
library(ggvenn)
library(ggplot2)

base_dir <- "~/Ovary_signatures/"

paths <- list(
  raw_full_deseq2  = file.path(base_dir, "0_DGE_raw",         "DE_full_OVARY_DESeq2.rds"),
  ae_full_deseq2   = file.path(base_dir, "1_DGE_autoencoder", "DE_full_OVARY_DESeq2.rds")
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

# ===================== FILTRADO (SIG + EFECTO) ===================== #
raw_sig <- DE_full_raw_DESeq2 %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF,
         !is.na(log2FoldChange), abs(log2FoldChange) > LFC_CUTOFF) %>%
  transmute(identifier, sign_raw = sign(log2FoldChange))

ae_sig <- DE_full_ae_DESeq2 %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF,
         !is.na(log2FoldChange), abs(log2FoldChange) > LFC_CUTOFF) %>%
  transmute(identifier, sign_ae = sign(log2FoldChange))

# ===================== COMPARTIDOS MISMA DIRECCIÓN (QC) ===================== #
shared_same_direction <- inner_join(raw_sig, ae_sig, by = "identifier") %>%
  filter(sign_raw == sign_ae)

cat("padj <", PADJ_CUTOFF, " & |log2FC| >", LFC_CUTOFF, "\n")
cat("Genes significativos+efecto (ovary_gtex_control):", n_distinct(raw_sig$identifier), "\n")
cat("Genes significativos+efecto (Reference_control):",  n_distinct(ae_sig$identifier), "\n")
cat("Genes compartidos (misma dirección):",             n_distinct(shared_same_direction$identifier), "\n")

# ===================== LISTAS PARA VENN ===================== #
venn_list <- list(
  ovary_gtex_control = unique(raw_sig$identifier),
  Reference_control  = unique(ae_sig$identifier)
)

# ===================== PLOT ===================== #
venn_plot <- ggvenn(
  venn_list,
  fill_color      = c("#2A9D8F", "#E76F51"),  # <-- colores nuevos
  fill_alpha      = 0.55,
  stroke_size     = 0.9,
  set_name_size   = 6,
  text_size       = 0,        # quita números
  show_percentage = FALSE     # quita porcentajes
) +
  labs(
    title = "Significant differentially expressed genes in ovarian cancer vs controls (DESeq2)",
    subtitle = paste0("padj < ", PADJ_CUTOFF, " and |log2FC| > ", LFC_CUTOFF)
  ) +
  theme_void(base_size = 13) +                 # quita ejes/grids/fondo tipo ggplot
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none"
  )

# ===================== GUARDAR PDF ===================== #
out_pdf <- file.path(base_dir, "venn_DESeq2_ovary_gtex_control_vs_Reference_control_padj0.01_log2FC1.pdf")
ggsave(out_pdf, plot = venn_plot, device = cairo_pdf, width = 7, height = 6, units = "in")

cat("PDF guardado en:", out_pdf, "\n")
