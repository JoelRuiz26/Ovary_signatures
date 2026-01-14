#!/usr/bin/env Rscript

#¿Los mismos genes cambian en la misma dirección y 
#son considerados igual de importantes en ambos datasets?
#a validación de consistencia de señal y prioridad

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ===================== CONFIG ===================== #
base_dir <- "~/Ovary_signatures"
paths <- list(
  ae_deseq2 = file.path(base_dir, "0_DGE_raw", "DE_full_OVARY_DESeq2.rds"),
  geo_rra   = file.path(base_dir, "2_DEG_GEO", "GEO_ovarian_cancer_RRA_results.rds")
)

out_dir <- file.path(base_dir, "4_AE_GEO_validation_rank_by_padj_and_direction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 1000
set.seed(1)

# ===================== HELPERS ===================== #
read_rds_safe <- function(p) {
  if (!file.exists(p)) stop("File not found: ", p)
  readRDS(p)
}

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  x
}

perm_p_cor <- function(x, y, method = c("spearman", "kendall"), n_perm = 1000) {
  method <- match.arg(method)
  obs <- suppressWarnings(cor(x, y, method = method))
  perm <- replicate(n_perm, suppressWarnings(cor(x, sample(y), method = method)))
  p <- (1 + sum(abs(perm) >= abs(obs))) / (n_perm + 1)
  list(obs = obs, p = p)
}

perm_p_dir_concordance <- function(ae_dir, geo_dir, n_perm = 1000) {
  obs <- mean(ae_dir == geo_dir)
  perm <- replicate(n_perm, mean(ae_dir == sample(geo_dir)))
  p <- (1 + sum(perm >= obs)) / (n_perm + 1)
  list(obs = obs, p = p)
}

# ===================== LOAD ===================== #
DE_ae <- read_rds_safe(paths$ae_deseq2)
GEO   <- read_rds_safe(paths$geo_rra)

# ===================== CHECK COLS ===================== #
req_ae  <- c("Symbol", "log2FoldChange", "padj")
req_geo <- c("Gene.symbol", "direction", "p_adj")

missing_ae  <- setdiff(req_ae,  colnames(DE_ae))
missing_geo <- setdiff(req_geo, colnames(GEO))

if (length(missing_ae)  > 0) stop("Missing columns in AE: ", paste(missing_ae, collapse = ", "))
if (length(missing_geo) > 0) stop("Missing columns in GEO: ", paste(missing_geo, collapse = ", "))

# ===================== STANDARDIZE + RANK (WITHIN DIRECTION BY PADJ) ===================== #
ae_df <- DE_ae %>%
  mutate(
    gene   = clean_gene(Symbol),
    log2FC = as.numeric(log2FoldChange),
    padj   = as.numeric(padj)
  ) %>%
  filter(!is.na(gene), !is.na(log2FC), !is.na(padj)) %>%
  mutate(direction_ae = ifelse(log2FC > 0, "up", "down")) %>%
  arrange(direction_ae, padj, desc(abs(log2FC)), gene) %>%
  distinct(gene, .keep_all = TRUE) %>%
  group_by(direction_ae) %>%
  mutate(rank_ae = row_number()) %>%
  ungroup()

geo_df <- GEO %>%
  mutate(
    gene      = clean_gene(Gene.symbol),
    direction = tolower(trimws(direction)),
    p_adj     = as.numeric(p_adj)
  ) %>%
  filter(!is.na(gene), !is.na(p_adj), direction %in% c("up", "down")) %>%
  arrange(direction, p_adj, gene) %>%
  distinct(gene, .keep_all = TRUE) %>%
  group_by(direction) %>%
  mutate(rank_geo = row_number()) %>%
  ungroup()

# ===================== MERGE ===================== #
m <- inner_join(
  ae_df  %>% select(gene, log2FC, padj, direction_ae, rank_ae),
  geo_df %>% select(gene, direction, p_adj, rank_geo),
  by = "gene"
)

if (nrow(m) < 30) warning("Few shared genes (n=", nrow(m), "). Interpret with caution.")

# ===================== DIRECTION AGREEMENT + MCC ===================== #
#(Matthew’s Correlation Coefficient)
#medida de asociación entre dos clasificaciones binarias.

#“¿Los genes más significativos en AE también lo son en GEO?”


m <- m %>% mutate(concordant_dir = (direction_ae == direction))
tab <- table(AE = m$direction_ae, GEO = m$direction)

# robust extraction (missing cells -> 0) + cast to numeric to avoid integer overflow
get_cell <- function(tt, r, c) {
  rn <- rownames(tt); cn <- colnames(tt)
  if (r %in% rn && c %in% cn) as.numeric(tt[r, c]) else 0
}
a <- get_cell(tab, "up",   "up")
b <- get_cell(tab, "up",   "down")
c <- get_cell(tab, "down", "up")
d <- get_cell(tab, "down", "down")

mcc_den <- sqrt((a + b) * (a + c) * (d + b) * (d + c))
MCC <- ifelse(mcc_den == 0, NA_real_, (a * d - b * c) / mcc_den)

dir_perm <- perm_p_dir_concordance(m$direction_ae, m$direction, N_PERM)

# ===================== RANK AGREEMENT (ONLY CONCORDANT DIRECTION; UP and DOWN SEPARATE) ===================== #
m_up   <- m %>% filter(direction_ae == "up",   direction == "up")
m_down <- m %>% filter(direction_ae == "down", direction == "down")

sp_up <- if (nrow(m_up) >= 30) perm_p_cor(m_up$rank_ae, m_up$rank_geo, "spearman", N_PERM) else list(obs = NA_real_, p = NA_real_)
kd_up <- if (nrow(m_up) >= 30) perm_p_cor(m_up$rank_ae, m_up$rank_geo, "kendall",  N_PERM) else list(obs = NA_real_, p = NA_real_)

sp_dn <- if (nrow(m_down) >= 30) perm_p_cor(m_down$rank_ae, m_down$rank_geo, "spearman", N_PERM) else list(obs = NA_real_, p = NA_real_)
kd_dn <- if (nrow(m_down) >= 30) perm_p_cor(m_down$rank_ae, m_down$rank_geo, "kendall",  N_PERM) else list(obs = NA_real_, p = NA_real_)

# ===================== PLOTS ===================== #
p_rank <- ggplot(
  bind_rows(
    m_up   %>% mutate(group = "UP (AE & GEO)"),
    m_down %>% mutate(group = "DOWN (AE & GEO)")
  ),
  aes(x = rank_ae, y = rank_geo)
) +
  geom_point(alpha = 0.35, size = 1.1) +
  scale_x_continuous(trans = "reverse") +
  scale_y_continuous(trans = "reverse") +
  facet_wrap(~ group, scales = "free") +
  labs(
    title = "Rank agreement within direction (ranked by adjusted p-value)",
    x = "AE rank (padj asc within direction)",
    y = "GEO rank (p_adj asc within direction)"
  ) +
  theme_classic(base_size = 14)

ggsave(file.path(out_dir, "4_0_rank_scatter_within_direction_AE_vs_GEO.pdf"),
       p_rank, width = 9.5, height = 5.8)

p_dir <- m %>%
  count(concordant_dir) %>%
  mutate(label = ifelse(concordant_dir, "Concordant", "Discordant")) %>%
  ggplot(aes(x = label, y = n)) +
  geom_col() +
  labs(
    title = "Direction agreement: AE vs GEO",
    x = NULL,
    y = "Number of shared genes"
  ) +
  theme_classic(base_size = 14)

ggsave(file.path(out_dir, "4_1_direction_counts.pdf"),
       p_dir, width = 7.5, height = 5.5)

# ===================== SUMMARY ===================== #
cat("Done.\n")
cat("Outputs: ", out_dir, "\n", sep = "")
cat("Shared genes: ", nrow(m), "\n", sep = "")

cat("Direction concordance rate: ",
    round(mean(m$concordant_dir), 4),
    " | perm p=", format(dir_perm$p, scientific = TRUE), "\n", sep = "")

cat("MCC: ", ifelse(is.na(MCC), "NA", round(MCC, 4)), "\n", sep = "")

cat("\nWithin-direction rank correlations (only genes concordant in direction):\n")

cat("UP genes: n=", nrow(m_up),
    " | Spearman=", round(sp_up$obs, 4), " p=", format(sp_up$p, scientific = TRUE),
    " | Kendall=",  round(kd_up$obs, 4), " p=", format(kd_up$p, scientific = TRUE), "\n", sep = "")

cat("DOWN genes: n=", nrow(m_down),
    " | Spearman=", round(sp_dn$obs, 4), " p=", format(sp_dn$p, scientific = TRUE),
    " | Kendall=",  round(kd_dn$obs, 4), " p=", format(kd_dn$p, scientific = TRUE), "\n", sep = "")

cat("\n2x2 table (AE_direction vs GEO_direction):\n")
print(tab)


#     GEO
#AE   down   up
#down 6089 1521
#up   3813 3085


#Spearman
#Correlación entre rangos
#Basada en diferencias cuadráticas de rangos
#“Si un gen está más arriba (más significativo) en AE, 
#¿tiende también a estar más arriba en GEO?”


#Kendall
#Medida de concordancia de pares
#“Si tomo dos genes al azar, ¿AE y GEO están de acuerdo en cuál es más significativo?”


#Done.
#Shared genes: 14508
#Direction concordance rate: 0.6177 | perm p=9.99001e-04
#MCC: 0.2481

#Within-direction rank correlations (only genes concordant in direction):
#  UP genes: n=2331 | Spearman=0.2636 p=9.99001e-04 | Kendall=0.1801 p=9.99001e-04
#DOWN genes: n=6631 | Spearman=0.1113 p=9.99001e-04 | Kendall=0.0743 p=9.99001e-04

#2x2 table (AE_direction vs GEO_direction):
#  GEO
#AE     down   up
#down  6631  988
#up    4558 2331

#Se compararon 14,508 genes comunes entre ambos datasets.
#La dirección del cambio (up/down) coincide en ≈62% de los genes, una proporción significativamente mayor a la esperada por azar (p ≈ 0.001 por permutaciones).
#El MCC = 0.25 indica una asociación positiva moderada entre las asignaciones de dirección en ambos datasets (no independencia, pero lejos de concordancia perfecta).

#Al evaluar solo genes con la misma dirección en ambos datasets:
#  Genes UP (n = 2,331) muestran una correlación positiva moderada en el ranking por significancia ajustada (Spearman ≈ 0.26; Kendall ≈ 0.18; p ≈ 0.001).
#Genes DOWN (n = 6,631) muestran una correlación positiva débil en el ranking (Spearman ≈ 0.11; Kendall ≈ 0.07; p ≈ 0.001), con mayor variabilidad entre datasets.
#Interpretación en una frase:
#  Existe una concordancia estadística no aleatoria pero parcial entre ambos datasets a nivel de genes: la dirección del cambio es consistentemente similar en una mayoría de genes, mientras que la priorización por significancia es más estable para genes upregulados que para genes downregulados, s








