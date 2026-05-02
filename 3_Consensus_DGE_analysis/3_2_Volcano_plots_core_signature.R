
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(vroom)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"

paths <- list(
  raw  = file.path(base_dir, "0_DGE_GTEx", "DE_full_OVARY_DESeq2_GTEx.rds"),
  ae   = file.path(base_dir, "1_DGE_AE",   "DE_full_OVARY_DESeq2_AE.rds"),
  core = "/STORAGE/csbig/jruiz/Ovary_data/3_Consensus_DGE_analysis/3_0_1_genes_concordant.tsv"
)

out_dir <- file.path(base_dir, "3_Consensus_DGE_analysis/3_2_DGE_volcano_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cut <- 0.01
lfc_cut  <- 1
TOP_N    <- 5   # <-- 5 Up y 5 Down por dataset

max_neglog10 <- 301
p_floor <- 1e-300

# ===================== HELPERS ===================== #
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("No existe: ", p)
  readRDS(p)
}

prep_df <- function(df) {
  
  if (!("Symbol" %in% colnames(df))) {
    if ("Symbol_autho" %in% colnames(df)) {
      df <- df %>% mutate(Symbol = as.character(Symbol_autho))
    } else {
      stop("No encuentro Symbol ni Symbol_autho.")
    }
  }
  
  padj_col <- if ("padj_safe" %in% colnames(df)) "padj_safe" else if ("padj" %in% colnames(df)) "padj" else NA
  if (is.na(padj_col)) stop("No encuentro padj_safe ni padj.")
  
  df %>%
    transmute(
      Gene = as.character(Symbol),
      log2FoldChange = as.numeric(log2FoldChange),
      padj_used = as.numeric(.data[[padj_col]])
    ) %>%
    filter(
      !is.na(Gene), Gene != "",
      is.finite(log2FoldChange),
      is.finite(padj_used),
      padj_used > 0
    ) %>%
    mutate(
      padj_used = pmax(padj_used, p_floor),
      neglog10  = -log10(padj_used)
    ) %>%
    filter(neglog10 <= max_neglog10) %>%
    mutate(
      Type = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Overexpressed",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Underexpressed",
        TRUE ~ "Not significant"
      ),
      direction = case_when(
        padj_used < padj_cut & log2FoldChange >=  lfc_cut ~ "Up",
        padj_used < padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ NA_character_
      ),
      score = abs(log2FoldChange) * neglog10
    )
}

# ---- pick Top genes to label: TOP_N Up + TOP_N Down (solo core_signature) ----
pick_top_labels_core <- function(df, core_tbl, top_n = 5) {
  
  cand <- df %>%
    inner_join(core_tbl, by = "Gene") %>%
    filter(!is.na(direction)) %>%
    filter(!grepl("-", Gene) & !grepl("\\.", Gene)) %>%
    group_by(Gene) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  top_up <- cand %>%
    filter(direction == "Up") %>%
    slice_max(order_by = score, n = top_n, with_ties = FALSE)
  
  top_down <- cand %>%
    filter(direction == "Down") %>%
    slice_max(order_by = score, n = top_n, with_ties = FALSE)
  
  bind_rows(top_up, top_down)
}

# ===================== RUN (SOLO GTEx + AE) ===================== #
raw_df <- prep_df(read_rds_safe(paths$raw)) %>% mutate(dataset = "GTEx control")
ae_df  <- prep_df(read_rds_safe(paths$ae))  %>% mutate(dataset = "AE control")

all_df <- bind_rows(raw_df, ae_df) %>%
  mutate(dataset = factor(dataset, levels = c("GTEx control", "AE control")))

# eje X común para ambos panels
x_max <- max(abs(all_df$log2FoldChange), na.rm = TRUE)
# eje X fijo para ambos panels (LFC -20 a 20)
x_lim <- c(-10, 10)

# margen extra para labels (evita que se recorten con nudge_x)
x_pad <- 4
x_lim_plot <- x_lim + c(-x_pad, x_pad)

# ---- core_signature (col gene) ----
if (!file.exists(paths$core)) stop("No existe core_signature: ", paths$core)
core_signature <- vroom(paths$core, show_col_types = FALSE)

if (!("gene" %in% colnames(core_signature))) {
  stop("core_signature no tiene columna 'gene'. Columnas: ", paste(colnames(core_signature), collapse = ", "))
}

core_tbl <- core_signature %>%
  filter(n_sources == 3) %>%
  transmute(Gene = as.character(gene)) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct()

all_df <- all_df %>%
  mutate(
    in_core = Gene %in% core_tbl$Gene,
    plot_group = case_when(
      in_core & Type == "Overexpressed"    ~ "Core_Up",
      in_core & Type == "Underexpressed"   ~ "Core_Down",
      !in_core & Type == "Not significant" ~ "NS",
      TRUE ~ "NonCore_Sig"
    )
  )

# ---- label df: top (core_signature) por dataset ----
label_df <- all_df %>%
  group_by(dataset) %>%
  group_modify(~ pick_top_labels_core(.x, core_tbl, top_n = TOP_N)) %>%
  ungroup()

# ===================== PLOT ===================== #
library(ggh4x)

p <- ggplot(all_df, aes(log2FoldChange, neglog10)) +
  geom_point(aes(color = plot_group, alpha = plot_group), size = 0.8) +
  scale_color_manual(
    values = c(
      "Core_Up"      = "#D7301F",
      "Core_Down"    = "#2C7FB8",
      "NS"           = "grey",
      "NonCore_Sig"  = "grey"
    ),
    breaks = c("Core_Up", "Core_Down"),
    labels = c("Overexpressed", "Underexpressed")
  ) +
  scale_alpha_manual(
    values = c(
      "Core_Up"     = 0.90,
      "Core_Down"   = 0.90,
      "NS"          = 0.60,
      "NonCore_Sig" = 0.1
    ),
    guide = "none"
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut),     linetype = "dashed", linewidth = 0.4) +
  geom_point(
    data = label_df,
    aes(x = log2FoldChange, y = neglog10, color = plot_group),
    inherit.aes = FALSE,
    size = 1.6,
    alpha = 0.95,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = label_df %>% filter(plot_group == "Core_Up"),
    aes(x = log2FoldChange, y = neglog10, label = Gene),
    inherit.aes = FALSE,
    size = 3.2,
    color = "#D7301F",
    segment.color = "#D7301F",
    segment.size = 0.30,
    box.padding = 0.20,
    point.padding = 0.15,
    force = 1.0,
    force_pull = 3.0,
    max.iter = 20000,
    max.time = 2,
    min.segment.length = 0.15,
    direction = "both",
    nudge_x = 0.35,
    seed = 132
  ) +
  ggrepel::geom_text_repel(
    data = label_df %>% filter(plot_group == "Core_Down"),
    aes(x = log2FoldChange, y = neglog10, label = Gene),
    inherit.aes = FALSE,
    size = 3.2,
    color = "#2C7FB8",
    segment.color = "#2C7FB8",
    segment.size = 0.30,
    box.padding = 0.20,
    point.padding = 0.15,
    force = 1.0,
    force_pull = 3.0,
    max.iter = 20000,
    max.time = 2,
    min.segment.length = 0.15,
    direction = "both",
    nudge_x = -0.35,
    seed = 32
  ) +
  ggh4x::facet_grid2(
    dataset ~ .,
    scales = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = c("#B5C3E3", "#F9EDB2"),
        colour = NA
      )
    )
  ) +
  coord_cartesian(xlim = x_lim_plot, ylim = c(0, max_neglog10 * 1.08), clip = "off") +
  labs(
    x = "log2FoldChange",
    y = "-log10(padj)",
    color = "Expression"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text.y = element_text(face = "bold", size = 12),
    strip.switch.pad.grid = unit(0.45, "cm"),
    panel.spacing.y = unit(0.35, "lines"),
    plot.margin = margin(10, 35, 10, 35)
  )

p

# ===================== SAVE ===================== #
out_pdf <- file.path(out_dir, "3_2_1_volcano_GTEx_AE_only_CORE.pdf")
out_png <- sub("\\.pdf$", ".png", out_pdf)

n_panels <- length(levels(all_df$dataset))
ggsave(out_pdf, p, width = 9.5, height = max(6, 2.3 * n_panels), units = "in", device = cairo_pdf)
ggsave(out_png, p, width = 9.5, height = max(6, 2.3 * n_panels), units = "in", dpi = 600, device = "png", bg = "white")

message("Done. Outputs en: ", out_dir)
