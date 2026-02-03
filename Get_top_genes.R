#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(writexl)
  library(httr)
  library(jsonlite)
  library(tibble)
})

options(stringsAsFactors = FALSE)

# -------------------- Paths -------------------- #
in_file  <- "/home/joelr/Ovary_signatures/3_Consensus_DGE_analysis/3_0_1_genes_concordant.tsv"
out_file <- "/home/joelr/Ovary_signatures/3_Consensus_DGE_analysis/3_0_1_genes_concordant_filtered.xlsx"

# -------------------- Read -------------------- #
df <- read_tsv(in_file, show_col_types = FALSE)

# -------------------- Genes vía adenosina (imagen) -------------------- #
genes_img <- c("ADA", "ADORA1", "ADORA2A", "ADORA3")

# -------------------- Filter + columna de interés -------------------- #
df_filtered <- df %>%
  filter(n_sources == 3 | gene %in% genes_img) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(Adenosine = if_else(gene %in% genes_img, "Adenosine", "other"))

# ============================================================
# ANOTACIÓN AUTOMÁTICA DESDE MyGene.info (sin Bioconductor)
# - canonical_name: nombre del gen
# - canonical_function: resumen/función (campo "summary" si existe)
# - reactome_pathways / kegg_pathways: si están disponibles
# ============================================================

collapse_pathways <- function(x) {
  # MyGene puede devolver: NULL / list / data.frame / vector
  if (is.null(x)) return(NA_character_)
  if (is.atomic(x)) return(paste(unique(as.character(x)), collapse = "; "))
  
  # Si es lista, intentamos extraer nombres comunes
  # Reactome suele venir como lista de objetos con "name"
  tryCatch({
    # Flatten a data.frame si se puede
    dfp <- as.data.frame(x)
    if ("name" %in% names(dfp)) {
      vals <- unique(as.character(dfp$name))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) == 0) return(NA_character_)
      return(paste(vals, collapse = "; "))
    }
    # si no hay name, colapsa todo a texto
    vals <- unique(unlist(x))
    vals <- vals[!is.na(vals) & vals != ""]
    if (length(vals) == 0) return(NA_character_)
    paste(as.character(vals), collapse = "; ")
  }, error = function(e) NA_character_)
}

query_mygene_batch <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) return(tibble())
  
  url <- "https://mygene.info/v3/query"
  fields <- "symbol,name,summary,pathway.reactome,pathway.kegg"
  
  # POST batch (mucho más eficiente que 5000 requests)
  resp <- httr::POST(
    url = url,
    encode = "json",
    body = list(
      q = symbols,
      scopes = "symbol",
      fields = fields,
      species = "human"
    ),
    httr::timeout(120)
  )
  
  if (httr::status_code(resp) != 200) {
    stop("MyGene.info respondió con status: ", httr::status_code(resp))
  }
  
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  res <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  
  # 'res' suele ser lista de resultados (uno por query)
  out <- lapply(res, function(r) {
    # Cuando no hay hit, puede traer 'notfound' o lista vacía
    if (is.null(r) || isTRUE(r$notfound)) {
      return(tibble(
        gene = NA_character_,
        canonical_name = NA_character_,
        canonical_function = NA_character_,
        reactome_pathways = NA_character_,
        kegg_pathways = NA_character_
      ))
    }
    
    tibble(
      gene = if (!is.null(r$query)) as.character(r$query) else NA_character_,
      canonical_name = if (!is.null(r$name)) as.character(r$name) else NA_character_,
      canonical_function = if (!is.null(r$summary)) as.character(r$summary) else NA_character_,
      reactome_pathways = collapse_pathways(r$pathway$reactome),
      kegg_pathways = collapse_pathways(r$pathway$kegg)
    )
  }) %>% dplyr::bind_rows()
  
  # algunas veces gene puede venir NA si cambió formato; entonces intenta usar symbol
  # pero nos quedamos con gene=query, que es tu símbolo original
  out %>%
    filter(!is.na(gene) & gene != "") %>%
    distinct(gene, .keep_all = TRUE)
}

annot <- tryCatch(
  query_mygene_batch(df_filtered$gene),
  error = function(e) {
    message("⚠️ No se pudo anotar con MyGene.info: ", conditionMessage(e))
    tibble(
      gene = unique(df_filtered$gene),
      canonical_name = NA_character_,
      canonical_function = NA_character_,
      reactome_pathways = NA_character_,
      kegg_pathways = NA_character_
    )
  }
)

df_filtered <- df_filtered %>%
  left_join(annot, by = "gene")

# -------------------- Checks -------------------- #
cat("Rows original:", nrow(df), "\n")
cat("Rows filtrado :", nrow(df_filtered), "\n\n")

cat("Adenosina (con anotación):\n")
print(
  df_filtered %>%
    filter(Adenosine == "via_adenosina") %>%
    dplyr::select(gene, n_sources, consensus_direction, Adenosine,
           canonical_name, canonical_function, reactome_pathways, kegg_pathways)
)

# -------------------- Write Excel -------------------- #
write_xlsx(df_filtered, out_file)
cat("\nExcel generado en:\n", out_file, "\n")
