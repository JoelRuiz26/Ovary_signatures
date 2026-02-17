
#load("~/Ovary_signatures/4_ORA_clusterProfiler/ORA_IMAGE.RData")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

# ============================================================
# DOTPLOT COMPACTO Y ELEGANTE (CORREGIDO)
# ============================================================

parse_ratio <- function(x) {
  sapply(x, function(z) {
    if (is.na(z) || !nzchar(z)) return(NA_real_)
    parts <- strsplit(z, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_real_)
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

TOP_N_PER_ONTOLOGY <- 12

plot_df <- ora_results %>%
  mutate(
    Direction = factor(
      Direction,
      levels = c("UP", "DOWN"),
      labels = c("Overexpressed", "Underexpressed")),
    Ontology  = factor(Ontology,  levels = c("BP", "MF", "CC")),
    GeneRatio_num = parse_ratio(GeneRatio),
    mlog10FDR = -log10(p.adjust),
    Description_wrapped = str_wrap(Description, width = 36)
  ) %>%
  filter(!is.na(p.adjust), !is.na(Count)) %>%
  group_by(Direction, Ontology) %>%
  slice_min(order_by = p.adjust, n = TOP_N_PER_ONTOLOGY, with_ties = FALSE) %>%
  ungroup()

# ordenar términos dentro de cada panel
plot_df <- plot_df %>%
  group_by(Ontology, Direction) %>%
  mutate(
    Description_wrapped = factor(
      Description_wrapped,
      levels = rev(unique(Description_wrapped[order(p.adjust)]))
    )
  ) %>%
  ungroup()

# ============================================================
# BUILD PLOT (aquí estaba el error)
# ============================================================

p_dot <- ggplot(plot_df, aes(x = Direction, y = Description_wrapped)) +
  geom_point(
    aes(size = Count, color = mlog10FDR),
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "D",
    direction = 1,
    breaks = seq(2, 10, by = 2),
    labels = seq(2, 10, by = 2),
    name = expression(-log[10]("FDR"))
  ) +
  scale_size(
    range = c(2.5, 8),
    name = "Gene count"
  ) +
  facet_grid(
    Ontology ~ .,
    scales = "free_y",
    space  = "free_y",
    switch = "y"
  ) +
  labs(
    title = "GO ORA (top terms; semantic-reduced)",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey95", color = "grey40"),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    panel.spacing.y = unit(0.6, "lines")
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(3.8, "cm"),
      barwidth  = unit(0.45, "cm"),
      ticks.colour = "grey20",
      frame.colour = "grey40"
    )
  )

p_dot <- p_dot +
  scale_color_viridis_c(
    option = "D",
    direction = 1,
    limits = c(2, 10),
    breaks = c(2, 4, 6, 8, 10),
    labels = c("2", "4", "6", "8", "10"),
    oob = scales::squish,          # <<< CLAVE
    name = expression(-log[10]~"(FDR)")
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(4.5, "cm"),
      barwidth  = unit(0.6, "cm"),
      ticks.colour = "grey20",
      frame.colour = "grey40"
    )
  )

# ============================================================
# SAVE
# ============================================================

ggsave(
  file.path(out_dir, "ORA_GO_topTerms_dotplot_COMPACT.png"),
  p_dot,
  width  = 7.8,
  height = 6.8,
  dpi    = 600
)

ggsave(
  file.path(out_dir, "ORA_GO_topTerms_dotplot_COMPACT.pdf"),
  p_dot,
  width  = 7.8,
  height = 6.8
)

cat("Compact dotplot saved.\n")
