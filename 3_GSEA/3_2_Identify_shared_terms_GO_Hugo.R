#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggvenn)
  library(stringr)
})

options(stringsAsFactors = FALSE)

# ==========================================================
# Script: Shared GO BP terms (GSEA outputs) — DESeq2 only
# - Uses existing TSV outputs in: 3_GSEA/3_GSEA_GO_HUGO
# - Compares:
#     * Ovary_control  = raw
#     * Reference_control = autoencoder
# - Shared terms: intersection by Description (no direction constraint)
# - Outputs saved in: 3_GSEA/3_GSEA_GO_HUGO/shared_DESeq2
#   * shared GO terms tables
#   * Venn diagram
#   * Dotplot Top20 Up+Down for each comparison (shared only)
# ==========================================================

# ===================== CONFIG ===================== #
go_dir <- "~/Ovary_signatures/3_GSEA/3_GSEA_GO_HUGO"
out_dir <- file.path(go_dir, "shared_DESeq2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define group names requested by user:
# - ovary_control = Raw
# - reference_control = Autoencoder
ovary_label <- "Raw"
ref_label   <- "Autoencoder"

# ===================== HELPERS ===================== #
pick_file <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No se encontró archivo con patrón: ", pattern, " en: ", dir)
  if (length(files) > 1) message("Más de un archivo coincide con patrón: ", pattern, "\nUsando: ", files[1])
  files[1]
}

load_gsea_go <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  needed <- c("ID", "Description", "NES", "p.adjust")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop("Faltan columnas requeridas en ", basename(path), ": ", paste(miss, collapse = ", "))
  df
}

get_direction <- function(df) {
  df %>% dplyr::mutate(Direction = ifelse(NES > 0, "Up", "Down"))
}

make_venn <- function(terms1, terms2, title_txt, subtitle_txt, out_pdf) {
  
  venn_list <- list(
    "Ovary_control"     = terms1,
    "Reference_control" = terms2
  )
  
  p <- ggvenn::ggvenn(
    venn_list,
    fill_color      = c("#1F77B4", "#FF7F0E"),
    fill_alpha      = 0.45,
    stroke_size     = 0.9,
    set_name_size   = 6,
    text_size       = 6,
    show_percentage = TRUE,
    digits          = 1,
    label_sep       = "\n"
  ) +
    ggplot2::labs(
      title = stringr::str_wrap(title_txt, 55),
      subtitle = subtitle_txt
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.line  = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  
  ggplot2::ggsave(out_pdf, p, width = 6, height = 6)
}

plot_top20_shared <- function(df, title, subtitle, out_pdf) {
  
  dfp <- df %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust), p.adjust > 0) %>%
    dplyr::mutate(neglog10_fdr = -log10(p.adjust)) %>%
    dplyr::arrange(NES)
  
  top_up <- dfp %>%
    dplyr::filter(NES > 0) %>%
    dplyr::arrange(p.adjust, dplyr::desc(NES)) %>%
    dplyr::slice_head(n = 20)
  
  top_down <- dfp %>%
    dplyr::filter(NES < 0) %>%
    dplyr::arrange(p.adjust, NES) %>%  # most negative first among downs
    dplyr::slice_head(n = 20)
  
  top <- dplyr::bind_rows(top_up, top_down) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot2::ggplot(top, ggplot2::aes(
    x = NES,
    y = Description,
    color = neglog10_fdr,
    size  = abs(NES)
  )) +
    ggplot2::geom_point(shape = 16, alpha = 0.95) +
    ggplot2::scale_color_viridis_c(name = expression(-log[10]("FDR"))) +
    ggplot2::scale_size_continuous(name = "|NES|") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  ggplot2::ggsave(out_pdf, p, width = 11, height = 8.5)
}

# ===================== INPUTS (DESeq2 ONLY) ===================== #
# Ovary_control = raw  (per your instruction)
raw_file <- pick_file(
  go_dir,
  pattern = "^GSEA_GO_BP_DGE_raw_DE_full_OVARY__DESeq2\\.tsv$"
)

# Reference_control = autoencoder (per your instruction)
ae_file <- pick_file(
  go_dir,
  pattern = "^GSEA_GO_BP_DGE_autoencoder_DE_full_OVARY__DESeq2\\.tsv$"
)

message("Using files:")
message("  Ovary_control (Raw): ", basename(raw_file))
message("  Reference_control (Autoencoder): ", basename(ae_file))

ovary <- load_gsea_go(raw_file) %>% get_direction()
ref   <- load_gsea_go(ae_file)  %>% get_direction()

# ===================== SHARED (NO direction constraint) ===================== #
shared_terms <- ovary %>%
  dplyr::inner_join(ref %>% dplyr::select(Description), by = "Description") %>%
  dplyr::distinct(Description)

# Save list of shared terms
readr::write_tsv(
  shared_terms,
  file.path(out_dir, "GO_BP_shared_terms_DESeq2.tsv")
)

# Save merged (lets you see both NES/p.adjust and direction for each side)
shared_merged <- ovary %>%
  dplyr::inner_join(ref, by = "Description", suffix = c(paste0("_", ovary_label), paste0("_", ref_label)))

readr::write_tsv(
  shared_merged,
  file.path(out_dir, "GO_BP_shared_merged_DESeq2.tsv")
)

# ===================== VENN ===================== #
make_venn(
  terms1 = ovary$Description,
  terms2 = ref$Description,
  title_txt = "GO Biological Process terms enriched in ovarian cancer",
  subtitle_txt = "Shared significant GO BP terms",
  out_pdf = file.path(out_dir, "GO_BP_shared_venn_DESeq2.pdf")
)

# ===================== DOTPLOTS (shared only; one per comparison) ===================== #
ovary_shared <- ovary %>%
  dplyr::semi_join(shared_terms, by = "Description")

ref_shared <- ref %>%
  dplyr::semi_join(shared_terms, by = "Description")

plot_top20_shared(
  df = ovary_shared,
  title = "Shared GO BP terms",
  subtitle = "Ovary tumors vs ovary control",
  out_pdf = file.path(out_dir, "GO_BP_shared_top20_up_down_Ovary_control_Raw_DESeq2.pdf")
)

plot_top20_shared(
  df = ref_shared,
  title = "Shared GO BP terms",
  subtitle = "Ovary tumors vs reference control",
  out_pdf = file.path(out_dir, "GO_BP_shared_top20_up_down_Reference_control_Autoencoder_DESeq2.pdf")
)

message("\nDone. Outputs in: ", out_dir)
