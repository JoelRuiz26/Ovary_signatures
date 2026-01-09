#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggvenn)
  library(stringr)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures/4_GSEA"

# Carpeta reales (según tu tree)
dbs <- list(
  KEGG         = "4_0_KEGG",
  Reactome     = "4_0_REACTOME",
  WikiPathways = "4_0_WIKIPATHWAYS"
)

# ===================== HELPERS ===================== #
load_sig <- function(dir, prefix) {
  
  # 1) Primero intenta con tu naming actual: *_SIG_0.01.tsv
  files <- list.files(
    path = dir,
    pattern = paste0("^", prefix, ".*_SIG_0\\.01\\.tsv$"),
    full.names = TRUE
  )
  
  # 2) Fallback por compatibilidad: *_0.01.tsv (ej. Reactome antiguo)
  if (length(files) == 0) {
    files <- list.files(
      path = dir,
      pattern = paste0("^", prefix, ".*_0\\.01\\.tsv$"),
      full.names = TRUE
    )
  }
  
  if (length(files) == 0) {
    stop("No se encontró archivo SIG para: ", prefix, " en ", dir)
  }
  if (length(files) > 1) {
    message("Más de un archivo encontrado para ", prefix, " en ", dir, "\nUsando: ", files[1])
  }
  
  readr::read_tsv(files[1], show_col_types = FALSE)
}

get_direction <- function(df) {
  df %>% dplyr::mutate(Direction = ifelse(NES > 0, "Up", "Down"))
}

clean_description <- function(df, db_name) {
  # Limpieza solicitada:
  # - Reactome: quitar "REACTOME_"
  # - WikiPathways: quitar "WP_"
  if (db_name == "Reactome") {
    df <- df %>% dplyr::mutate(Description = gsub("^REACTOME_", "", Description))
  }
  if (db_name == "WikiPathways") {
    df <- df %>% dplyr::mutate(Description = gsub("^WP_", "", Description))
  }
  df
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
  
  # PDF
  ggplot2::ggsave(out_pdf, p, width = 6, height = 6)
  
  # PNG 600 DPI
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggplot2::ggsave(out_png, p, width = 6, height = 6, dpi = 600)
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
    dplyr::arrange(p.adjust, NES) %>%
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
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold")
    )
  
  # PDF
  ggplot2::ggsave(out_pdf, p, width = 11, height = 8.5)
  
  # PNG 600 DPI
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggplot2::ggsave(out_png, p, width = 11, height = 8.5, dpi = 600)
}

# ===================== MAIN LOOP ===================== #
for (db in names(dbs)) {
  
  message("Processing ", db, "...")
  
  db_dir <- file.path(base_dir, dbs[[db]])
  shared_dir <- file.path(db_dir, "shared")
  dir.create(shared_dir, showWarnings = FALSE)
  
  ovary <- load_sig(db_dir, "GTEx") %>%
    clean_description(db) %>%   # <- NO: esto pasa db como df
    get_direction()
  
  # CORRECCIÓN: clean_description(., db)
  ovary <- load_sig(db_dir, "GTEx") %>%
    clean_description(., db) %>%
    get_direction()
  
  ref <- load_sig(db_dir, "AE") %>%
    clean_description(., db) %>%
    get_direction()
  
  # ---- compartidos (SIN exigir misma dirección) ----
  shared_terms <- ovary %>%
    dplyr::inner_join(
      ref %>% dplyr::select(Description),
      by = "Description"
    ) %>%
    dplyr::distinct(Description)
  
  # guardar lista de términos compartidos
  readr::write_tsv(
    shared_terms,
    file.path(shared_dir, paste0(db, "_shared_terms.tsv"))
  )
  
  # merged completo (útil para inspección; aquí sí verás dirección y stats de ambos)
  shared_merged <- ovary %>%
    dplyr::inner_join(ref, by = "Description", suffix = c("_GTEx", "_AE"))
  
  readr::write_tsv(
    shared_merged,
    file.path(shared_dir, paste0(db, "_shared_merged.tsv"))
  )
  
  # ---- Venn (compartidos, sin dirección) ----
  make_venn(
    terms1 = ovary$Description,
    terms2 = ref$Description,
    title_txt = paste(db, "pathways enriched in ovarian cancer"),
    subtitle_txt = "Shared significant pathways",
    out_pdf = file.path(shared_dir, paste0(db, "_shared_venn.pdf"))
  )
  
  # ---- Dotplot 1: GTEx (solo compartidos) ----
  ovary_shared <- ovary %>%
    dplyr::semi_join(shared_terms, by = "Description")
  
  plot_top20_shared(
    df = ovary_shared,
    title = paste("Shared", db, "pathways"),
    subtitle = "Ovary tumors vs ovary control GTEx",
    out_pdf = file.path(shared_dir, paste0(db, "_shared_top20_up_down_GTEx.pdf"))
  )
  
  # ---- Dotplot 2: AE (solo compartidos) ----
  ref_shared <- ref %>%
    dplyr::semi_join(shared_terms, by = "Description")
  
  plot_top20_shared(
    df = ref_shared,
    title = paste("Shared", db, "pathways"),
    subtitle = "Ovary tumors vs reference control",
    out_pdf = file.path(shared_dir, paste0(db, "_shared_top20_up_down_AE.pdf"))
  )
}

message("\nAll shared analyses completed.")
