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

# Carpetas reales (según tu tree)
dbs <- list(
  KEGG         = "4_0_KEGG",
  Reactome     = "4_0_REACTOME",
  WikiPathways = "4_0_WIKIPATHWAYS"
)

# Umbral actual en tus outputs
fdr_tag <- "0.05"

# ===================== HELPERS ===================== #
load_sig <- function(dir, prefix, fdr_tag = "0.05") {
  
  # 1) Naming actual: *_SIG_0.05.tsv
  files <- list.files(
    path = dir,
    pattern = paste0("^", prefix, ".*_SIG_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
    full.names = TRUE
  )
  
  # 2) Fallback por compatibilidad: *_0.05.tsv (naming viejo)
  if (length(files) == 0) {
    files <- list.files(
      path = dir,
      pattern = paste0("^", prefix, ".*_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
      full.names = TRUE
    )
  }
  
  if (length(files) == 0) {
    stop("No se encontró archivo SIG para: ", prefix, " (FDR ", fdr_tag, ") en ", dir)
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
    "Set_1" = terms1,
    "Set_2" = terms2
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
  
  ggplot2::ggsave(out_pdf, p, width = 11, height = 8.5)
  
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  ggplot2::ggsave(out_png, p, width = 11, height = 8.5, dpi = 600)
}

# ---- Core: genera "exactamente lo mismo" para un par A vs B ----
pairwise_shared <- function(db, shared_dir, A_prefix, A_label, A_subtitle,
                            B_prefix, B_label, B_subtitle, fdr_tag = "0.05") {
  
  # Cargar y limpiar
  A <- load_sig(file.path(base_dir, dbs[[db]]), A_prefix, fdr_tag = fdr_tag) %>%
    clean_description(., db) %>%
    get_direction()
  
  B <- load_sig(file.path(base_dir, dbs[[db]]), B_prefix, fdr_tag = fdr_tag) %>%
    clean_description(., db) %>%
    get_direction()
  
  # Compartidos (sin exigir misma dirección)
  shared_terms <- A %>%
    dplyr::inner_join(B %>% dplyr::select(Description), by = "Description") %>%
    dplyr::distinct(Description)
  
  # Guardar lista de términos compartidos
  readr::write_tsv(
    shared_terms,
    file.path(shared_dir, paste0(db, "_shared_terms_", A_prefix, "_vs_", B_prefix, ".tsv"))
  )
  
  # Merged completo
  shared_merged <- A %>%
    dplyr::inner_join(B, by = "Description",
                      suffix = c(paste0("_", A_label), paste0("_", B_label)))
  
  readr::write_tsv(
    shared_merged,
    file.path(shared_dir, paste0(db, "_shared_merged_", A_prefix, "_vs_", B_prefix, ".tsv"))
  )
  
  # Venn
  make_venn(
    terms1 = A$Description,
    terms2 = B$Description,
    title_txt = paste(db, "pathways enriched in ovarian cancer"),
    subtitle_txt = paste0("Shared significant pathways: ", A_label, " vs ", B_label),
    out_pdf = file.path(shared_dir, paste0(db, "_shared_venn_", A_prefix, "_vs_", B_prefix, ".pdf"))
  )
  
  # Dotplot A (solo compartidos)
  A_shared <- A %>% dplyr::semi_join(shared_terms, by = "Description")
  plot_top20_shared(
    df = A_shared,
    title = paste("Shared", db, "pathways"),
    subtitle = A_subtitle,
    out_pdf = file.path(shared_dir, paste0(db, "_shared_top20_up_down_", A_prefix, "_vs_", B_prefix, "_", A_label, ".pdf"))
  )
  
  # Dotplot B (solo compartidos)
  B_shared <- B %>% dplyr::semi_join(shared_terms, by = "Description")
  plot_top20_shared(
    df = B_shared,
    title = paste("Shared", db, "pathways"),
    subtitle = B_subtitle,
    out_pdf = file.path(shared_dir, paste0(db, "_shared_top20_up_down_", A_prefix, "_vs_", B_prefix, "_", B_label, ".pdf"))
  )
}

# ===================== MAIN LOOP ===================== #
for (db in names(dbs)) {
  
  message("Processing ", db, "...")
  
  db_dir <- file.path(base_dir, dbs[[db]])
  shared_dir <- file.path(db_dir, "shared")
  dir.create(shared_dir, showWarnings = FALSE)
  
  # ---------- 1) ORIGINAL: GTEx vs AE ----------
  pairwise_shared(
    db = db,
    shared_dir = shared_dir,
    A_prefix = "GTEx",
    A_label = "GTEx",
    A_subtitle = "Ovary tumors vs ovary control GTEx",
    B_prefix = "AE",
    B_label = "AE",
    B_subtitle = "Ovary tumors vs reference control",
    fdr_tag = fdr_tag
  )
  
  # ---------- 2) GEO outputs ----------
  geo_sig_files <- list.files(
    path = db_dir,
    pattern = paste0("^GEO_.*_SIG_", gsub("\\.", "\\\\.", fdr_tag), "\\.tsv$"),
    full.names = FALSE
  )
  
  if (length(geo_sig_files) == 0) {
    message("  No GEO SIG files found for ", db, " (FDR ", fdr_tag, "). Skipping GEO.")
    next
  }
  
  # ====== FIX MÍNIMO (CRÍTICO) ======
  # Antes: dependía de 'db' ("Reactome") vs filename ("_REACTOME_") -> no matcheaba
  # Ahora: extrae el prefijo "GEO_GSE####" robustamente para KEGG/REACTOME/WIKIPATHWAYS
  geo_prefixes <- unique(sub("_.*$", "", geo_sig_files))
  # Ejemplo: "GEO_GSE14407"
  
  for (geo_prefix in geo_prefixes) {
    message("  GEO: ", geo_prefix)
    
    # (a) GEO vs GTEx
    pairwise_shared(
      db = db,
      shared_dir = shared_dir,
      A_prefix = geo_prefix,
      A_label = "GEO",
      A_subtitle = paste0("Ovary tumors vs control (", sub("^GEO_", "", geo_prefix), ")"),
      B_prefix = "GTEx",
      B_label = "GTEx",
      B_subtitle = "Ovary tumors vs ovary control GTEx",
      fdr_tag = fdr_tag
    )
    
    # (b) GEO vs AE
    pairwise_shared(
      db = db,
      shared_dir = shared_dir,
      A_prefix = geo_prefix,
      A_label = "GEO",
      A_subtitle = paste0("Ovary tumors vs control (", sub("^GEO_", "", geo_prefix), ")"),
      B_prefix = "AE",
      B_label = "AE",
      B_subtitle = "Ovary tumors vs reference control",
      fdr_tag = fdr_tag
    )
  }
}

message("\nAll shared analyses completed (including GEO).")
