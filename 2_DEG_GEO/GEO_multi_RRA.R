setwd("~/Ovary_signatures/2_DEG_GEO")

## ============================================================
## Required packages
## ============================================================

library(GEOquery)
library(limma)
library(gprofiler2)
library(RobustRankAggreg)
library(dplyr)



## ============================================================
## 1) Descargar datos GEO: lista de ExpressionSet
## ============================================================
# gse_ids: vector of IDs, e.g. c("GSE14407","GSE18520",...)
# destdir: folder where GEO files are stored
download_gse_series <- function(gse_ids,
                                destdir  = "geo_data",
                                platform = "GPL570",
                                retries  = 2,
                                gse_matrix = TRUE) {
  if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
  
  eset_list <- list()
  
  for (gse in gse_ids) {
    message("\n=== Procesando ", gse, " ===")
    
    # Expected path of the series_matrix file
    sm_file <- file.path(destdir, paste0(gse, "_series_matrix.txt.gz"))
    
    attempt <- 1
    eset_ok <- NULL
    
    while (attempt <= retries && is.null(eset_ok)) {
      message("   Intento ", attempt, " de ", retries)
      
      # Si el archivo existe, GEOquery debería reutilizarlo sin bajar de nuevo.
      # Si está corrupto, getGEO fallará y lo manejamos en el catch.
      eset_obj <- tryCatch(
        {
          GEOquery::getGEO(gse,
                           GSEMatrix = gse_matrix,
                           getGPL   = TRUE,
                           destdir  = destdir)
        },
        error = function(e) {
          warning("   ERROR al obtener ", gse, " (intento ", attempt, "): ",
                  conditionMessage(e))
          # Delete potentially corrupt series_matrix file
          if (file.exists(sm_file)) {
            message("   Borrando posible archivo parcial: ", sm_file)
            try(unlink(sm_file), silent = TRUE)
          }
          NULL
        }
      )
      
      if (is.null(eset_obj)) {
        attempt <- attempt + 1
        next
      }
      
      ## Platform selection
      
      if (methods::is(eset_obj, "ExpressionSet")) {
        # Solo un ExpressionSet
        eset_candidates <- list(eset_obj)
      } else if (is.list(eset_obj)) {
        # getGEO devolvió lista de ExpressionSet (varias plataformas)
        eset_candidates <- eset_obj
      } else {
        warning("   Objeto devuelto para ", gse,
                " no es ExpressionSet ni lista. Lo omito.")
        break
      }
      
      annots <- vapply(eset_candidates, annotation, character(1))
      message("   Plataformas encontradas: ", paste(unique(annots), collapse = ", "))
      
      idx <- which(annots == platform)
      if (length(idx) == 0) {
        warning("   Ningún ExpressionSet de ", gse,
                " usa plataforma ", platform,
                ". Tomo el primero por default.")
        idx <- 1L
      } else if (length(idx) > 1) {
        warning("   Más de un ExpressionSet con plataforma ", platform,
                " en ", gse, ". Tomo el primero.")
        idx <- idx[1]
      }
      
      eset <- eset_candidates[[idx]]
      
      if (!methods::is(eset, "ExpressionSet")) {
        warning("   El objeto elegido para ", gse,
                " no es un ExpressionSet válido. Intento de nuevo.")
        attempt <- attempt + 1
        next
      }
      
      message("   Usando plataforma ", annotation(eset), " para ", gse, ".")
      eset_ok <- eset
    } # while
    
    if (is.null(eset_ok)) {
      warning("   No se pudo obtener un ExpressionSet válido para ", gse,
              " tras ", retries, " intentos. Lo omito.")
      next
    }
    
    eset_list[[gse]] <- eset_ok
  }
  
  return(eset_list)
}

## ============================================================
## 2) Explore “conditions” in phenoData
## ============================================================
# eset_list: named list of ExpressionSet objects (output of download_gse_series)
# max_levels: maximum number of distinct levels to consider a column a “condition” type
inspect_gse_conditions <- function(eset_list, max_levels = 12) {
  stopifnot(is.list(eset_list))
  
  summary_list <- list()
  
  for (gse in names(eset_list)) {
    eset <- eset_list[[gse]]
    pheno <- Biobase::pData(eset)
    
    message("\n===================================================")
    message("GSE: ", gse)
    message("N muestras: ", nrow(pheno))
    message("Columnas en phenoData:")
    print(colnames(pheno))
    
    # candidate columns for "condition"
    cand_cols <- colnames(pheno)[vapply(pheno, function(x) {
      # treat as informative factor
      ux <- unique(as.character(x))
      n_ux <- length(ux)
      n_ux > 1 && n_ux <= max_levels
    }, logical(1))]
    
    if (length(cand_cols) == 0) {
      message("No se encontraron columnas con ≤ ", max_levels,
              " niveles distintos (posibles condiciones).")
      summary_list[[gse]] <- NULL
      next
    }
    
    message("\nColumnas candidatas a 'condición' (≤ ", max_levels, " niveles):")
    print(cand_cols)
    
    # frequency tables per candidate column
    cond_summary <- lapply(cand_cols, function(col) {
      tbl <- table(pheno[[col]], useNA = "ifany")
      sort(tbl, decreasing = TRUE)
    })
    names(cond_summary) <- cand_cols
    
    # readable output
    for (col in cand_cols) {
      cat("\n---", gse, "::", col, "---\n")
      print(cond_summary[[col]])
    }
    
    summary_list[[gse]] <- cond_summary
  }
  
  invisible(summary_list)
}

## ============================================================
## 3) Validate condition specifications
## ============================================================

check_contrast_specs <- function(esets, contrast_specs) {
  out <- list()
  
  for (gse in names(contrast_specs)) {
    message("=== Revisando ", gse, " ===")
    
    if (!gse %in% names(esets)) {
      warning("  -> ", gse, " NO está en 'esets'")
      next
    }
    
    eset <- esets[[gse]]
    cs   <- contrast_specs[[gse]]
    
    pd <- Biobase::pData(eset)
    
    # 1) Check column
    if (!cs$group_col %in% colnames(pd)) {
      warning("  -> Columna '", cs$group_col, "' NO existe en pData(", gse, ")")
      next
    }
    
    groups <- pd[[cs$group_col]]
    
    # 2) Table of raw levels
    tbl <- table(groups, useNA = "ifany")
    message("  Niveles encontrados:")
    print(tbl)
    
    control_levels <- cs$control
    case_spec      <- cs$case
    
    # ---------------------------------------------------
    # 3) Build control and case masks
    # ---------------------------------------------------
    control_mask <- rep(FALSE, length(groups))

    # 3a) Explicit controls (non-NA)
    non_na_ctrl <- control_levels[!is.na(control_levels)]
    if (length(non_na_ctrl) > 0) {
      missing_ctrl <- setdiff(non_na_ctrl, names(tbl))
      if (length(missing_ctrl) > 0) {
        warning("  -> Controles NO encontrados: ",
                paste(missing_ctrl, collapse = ", "))
      }
      control_mask <- control_mask | (!is.na(groups) & groups %in% non_na_ctrl)
    }
    
    # 3b) Controls marked as NA in pData
    if (any(is.na(control_levels))) {
      control_mask <- control_mask | is.na(groups)
    }
    
    # 3c) Cases
    if (identical(case_spec, "leftover")) {
      # everything not assigned as control is a case
      case_mask <- !control_mask

      # For reporting: levels used as case (ignoring NA)
      case_levels_report <- setdiff(names(tbl), non_na_ctrl)
      # if non_na_ctrl is empty and control is only NA, take all except "<NA>"
      if (any(is.na(control_levels))) {
        case_levels_report <- setdiff(names(tbl), NA_character_)
      }
      
    } else {
      # casos específicos por etiqueta
      missing_case <- setdiff(case_spec, names(tbl))
      if (length(missing_case) > 0) {
        warning("  -> Casos NO encontrados: ",
                paste(missing_case, collapse = ", "))
      }
      case_mask <- (!is.na(groups) & groups %in% case_spec)
      case_levels_report <- case_spec
    }
    
    n_control <- sum(control_mask)
    n_case    <- sum(case_mask)
    
    # ---------------------------------------------------
    # 4) Readable messages
    # ---------------------------------------------------
    message("  -> Control levels declarados: ",
            paste(ifelse(is.na(control_levels), "NA", control_levels),
                  collapse = ", "))
    message("     n_control = ", n_control)
    
    if (identical(case_spec, "leftover")) {
      message("  -> Case = leftover")
      message("     Niveles usados como caso: ",
              paste(case_levels_report, collapse = ", "))
    } else {
      message("  -> Case levels declarados: ",
              paste(case_levels_report, collapse = ", "))
    }
    message("     n_case    = ", n_case, "\n")
    
    # ---------------------------------------------------
    # 5) Save diagnostics
    # ---------------------------------------------------
    out[[gse]] <- list(
      found_levels   = tbl,
      control_levels = control_levels,
      case_levels    = case_levels_report,
      n_control      = n_control,
      n_case         = n_case
    )
  }
  
  return(out)
}


## ============================================================
## 4) Differential expression with limma for each dataset
## ============================================================
# eset_list: list of ExpressionSet objects
# contrast_specs: named list per GSE, for example:
#   list(
#     GSE14407 = list(
#       group_col = "source_name_ch1",
#       case      = c("serous ovarian cancer epithelium"),
#       control   = c("normal ovarian surface epithelium")
#     ),
#     GSE18520 = list(
#       group_col = "characteristics_ch1",
#       case      = c("advanced stage, high-grade primary OC"),
#       control   = c("normal ovarian surface epithelium")
#     )
#   )
#
# gene_annot_cols: columns from fData(eset) to join to the DE table.
#                  typically c("ID","Gene.symbol","Gene.title")
# log2_auto: if TRUE, tries to determine whether log2 transformation is needed.
limma_de_all <- function(eset_list,
                         contrast_specs,
                         gene_annot_cols = c("ID", "Gene.symbol", "Gene.title"),
                         log2_auto = TRUE) {
  stopifnot(is.list(eset_list))
  stopifnot(is.list(contrast_specs))
  
  de_list <- list()
  
  for (gse in names(eset_list)) {
    message("\n=== limma para ", gse, " ===")
    
    if (is.null(contrast_specs[[gse]])) {
      warning("   No hay spec de contraste para ", gse,
              ". Se omite este dataset.")
      next
    }
    
    spec <- contrast_specs[[gse]]
    required <- c("group_col", "case", "control")
    if (!all(required %in% names(spec))) {
      warning("   contrast_specs[['", gse, "']] debe tener: ",
              paste(required, collapse = ", "), ". Se omite.")
      next
    }
    
    eset  <- eset_list[[gse]]
    expr  <- Biobase::exprs(eset)
    pheno <- Biobase::pData(eset)
    
    # =========================
    # 3.1 Build group vector
    # =========================
    if (!spec$group_col %in% colnames(pheno)) {
      warning("   Columna '", spec$group_col, "' no está en phenoData de ",
              gse, ". Se omite.")
      next
    }
    col_vals <- pheno[[spec$group_col]]
    
    message("   Resumen de '", spec$group_col, "' en ", gse, ":")
    print(table(col_vals, useNA = "ifany"))
    
    # Leftover logic consistent with check_contrast_specs
    if (length(spec$case) == 1L &&
        spec$case[1] %in% c("leftover", "others", "not_control")) {

      message("   Using 'case = leftover' mode: everything not assigned as control is a case.")

      ctrl_vals     <- spec$control
      non_na_ctrl   <- ctrl_vals[!is.na(ctrl_vals)]
      control_mask  <- rep(FALSE, length(col_vals))

      # 1) explicit controls (non-NA)
      if (length(non_na_ctrl) > 0) {
        control_mask <- control_mask |
          (!is.na(col_vals) & col_vals %in% non_na_ctrl)
      }

      # 2) controls defined as NA in pData
      if (any(is.na(ctrl_vals))) {
        control_mask <- control_mask | is.na(col_vals)
      }
      
      group <- ifelse(control_mask, "control", "case")
      
    } else {
      # Standard mode: explicit case/control labels
      col_chr <- as.character(col_vals)
      group <- ifelse(col_chr %in% spec$case, "case",
               ifelse(col_chr %in% spec$control, "control", NA))
    }
    
    if (any(is.na(group))) {
      warning("   Algunas muestras de ", gse,
              " no se asignaron a case/control (NA en group). Las excluyo.")
    }
    
    keep <- !is.na(group)
    expr  <- expr[, keep, drop = FALSE]
    group <- factor(group[keep], levels = c("control", "case"))
    
    message("   Tabla final de grupos (tras filtrar NA):")
    print(table(group))
    
    if (length(unique(group)) < 2L) {
      warning("   Solo hay un nivel en group para ", gse,
              " después del filtrado. Se omite.")
      next
    }
    
    # =========================
    # 3.2 Automatic log2 detection
    # =========================
    if (log2_auto) {
      rng <- range(expr, na.rm = TRUE)
      if (rng[2] > 100 || (rng[2] - rng[1]) > 50) {
        message("   Asumo datos en escala lineal. Aplico log2(expr + 1).")
        expr <- log2(expr + 1)
      } else {
        message("   Asumo datos ya en log2 (microarreglos tipo GPL570).")
      }
    }
    
    # =========================
    # 3.3 Design matrix and limma fit
    # =========================
    design <- model.matrix(~ group)  # intercepto + efecto case vs control
    fit    <- limma::lmFit(expr, design)
    fit    <- limma::eBayes(fit)
    
    tt <- limma::topTable(fit,
                          coef          = "groupcase",
                          number        = Inf,
                          adjust.method = "BH",
                          sort.by       = "P")
    
    # =========================
    # 3.4 Add probe annotations
    # =========================
    fdat <- Biobase::fData(eset)
    common_probes <- intersect(rownames(tt), rownames(fdat))
    
    if (length(common_probes) > 0) {
      annot <- fdat[common_probes,
                    gene_annot_cols[gene_annot_cols %in% colnames(fdat)],
                    drop = FALSE]
      tt <- cbind(annot[rownames(tt), , drop = FALSE], tt)
    } else {
      message("   No hubo intersección entre rownames(expr) y rownames(fData). ",
              "Devuelvo solo la tabla limma.")
    }
    
    de_list[[gse]] <- tt
  }
  
  return(de_list)
}

## ============================================================
## 5) Annotate
## ============================================================
add_gene_symbol_with_gconvert_safe <- function(de_list,
                                               id_col   = "ID",
                                               organism = "hsapiens") {
  stopifnot(is.list(de_list))
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Instala 'gprofiler2' primero.")
  }
  
  out <- vector("list", length(de_list))
  names(out) <- names(de_list)
  
  for (nm in names(de_list)) {
    message("Mapeando probes → símbolos para ", nm, " ...")
    tt <- de_list[[nm]]
    
    if (is.null(tt) || nrow(tt) == 0L) {
      warning("  -> ", nm, " vacío o NULL. Lo devuelvo tal cual.")
      out[[nm]] <- tt
      next
    }
    
    ids <- if (id_col %in% colnames(tt)) {
      as.character(tt[[id_col]])
    } else {
      rownames(tt)
    }
    
    # gconvert returns one row per query (NA in 'name' if no mapping)
    gc <- gprofiler2::gconvert(
      query      = ids,
      organism   = organism,
      target     = "HGNC",
      mthreshold = 1,
      filter_na  = FALSE
    )

    # Identify the query index column and the symbol column
    idx_col  <- colnames(gc)[1]          # usually "query_number" or similar
    sym_col  <- "name"                   # standard column with symbol / name

    # Reorder by query index (1..n)
    gc[[idx_col]] <- as.integer(gc[[idx_col]])
    gc <- gc[order(gc[[idx_col]]), , drop = FALSE]

    # Symbol vector aligned to 'ids' by position
    gene_symbol <- gc[[sym_col]]
    
    if (length(gene_symbol) != nrow(tt)) {
      warning("  -> Longitud de símbolos (", length(gene_symbol),
              ") != nrow(tt) (", nrow(tt), ") en ", nm,
              ". Algo raro pasó con gconvert.")
    } else {
      tt$Gene.symbol <- gene_symbol
    }
    
    out[[nm]] <- tt
  }
  
  out
}

add_gene_symbol_with_gconvert <- function(de_list,
                                          id_col    = "ID",
                                          organism  = "hsapiens",
                                          chunk_size = 1000L) {
  stopifnot(is.list(de_list))
  
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Necesito el paquete 'gprofiler2'. Instálalo con install.packages('gprofiler2').")
  }
  
  out <- vector("list", length(de_list))
  names(out) <- names(de_list)
  
  for (nm in names(de_list)) {
    message("Mapeando probes → símbolos con gconvert para ", nm, " ...")
    tt <- de_list[[nm]]
    
    if (is.null(tt) || nrow(tt) == 0L) {
      warning("  -> ", nm, " está vacío o es NULL. Lo devuelvo tal cual.")
      out[[nm]] <- tt
      next
    }
    
    # 1) Tomar IDs de probe
    if (id_col %in% colnames(tt)) {
      probe_ids <- as.character(tt[[id_col]])
    } else {
      message("  -> No encontré columna '", id_col,
              "'. Uso rownames(tt) como IDs.")
      probe_ids <- rownames(tt)
    }
    
    probe_ids <- probe_ids[!is.na(probe_ids) & nzchar(probe_ids)]
    unique_ids <- unique(probe_ids)
    
    if (!length(unique_ids)) {
      warning("  -> No hay IDs válidos en ", nm, ". Lo devuelvo tal cual.")
      out[[nm]] <- tt
      next
    }
    
    # 2) Split into chunks to avoid overloading the g:Profiler API
    idx <- seq_along(unique_ids)
    split_ids <- split(unique_ids, ceiling(idx / chunk_size))
    
    conv_list <- lapply(split_ids, function(v) {
      gprofiler2::gconvert(
        query      = v,
        organism   = organism,
        target     = "ENSG",  # el target da igual; usamos la columna 'name'
        mthreshold = Inf,
        filter_na  = TRUE
      )
    })
    
    conv <- do.call(rbind, conv_list)
    
    if (is.null(conv) || nrow(conv) == 0L) {
      warning("  -> gconvert no devolvió nada útil para ", nm,
              ". Devuelvo la tabla sin Gene.symbol.")
      out[[nm]] <- tt
      next
    }
    
    # 3) Expect columns 'input' (original ID) and 'name' (symbol)
    if (!all(c("input", "name") %in% colnames(conv))) {
      warning("  -> La salida de gconvert no tiene columnas 'input' y 'name' en ", nm,
              ". Revisa versión de gprofiler2.")
      out[[nm]] <- tt
      next
    }
    
    conv <- conv[!is.na(conv$input) & !is.na(conv$name), , drop = FALSE]
    conv <- conv[!duplicated(conv$input), , drop = FALSE]  # one symbol per input
    
    sym_map <- conv$name
    names(sym_map) <- conv$input  # input = ID de probe
    
    # 4) Map back to the order of tt
    ids_for_tt <- if (id_col %in% colnames(tt)) {
      as.character(tt[[id_col]])
    } else {
      rownames(tt)
    }
    
    gene_symbol <- sym_map[ids_for_tt]
    tt$Gene.symbol <- gene_symbol
    
    out[[nm]] <- tt
  }
  
  out
}

add_gene_symbol_from_fdata <- function(de_list,
                                       eset_list,
                                       out_col = "Gene.symbol",
                                       symbol_candidates = c(
                                         "Gene.symbol",
                                         "Gene Symbol",
                                         "GENE_SYMBOL",
                                         "Symbol",
                                         "SYMBOL"
                                       )) {
  stopifnot(is.list(de_list))
  stopifnot(is.list(eset_list))
  
  out <- list()
  
  for (nm in names(de_list)) {
    message("Añadiendo símbolos para ", nm, " ...")
    
    tt   <- de_list[[nm]]
    eset <- eset_list[[nm]]
    
    if (is.null(tt) || nrow(tt) == 0L) {
      warning("  -> Tabla vacía o NULL para ", nm, ". La devuelvo tal cual.")
      out[[nm]] <- tt
      next
    }
    if (is.null(eset)) {
      warning("  -> No hay ExpressionSet para ", nm, ". No puedo mapear símbolos.")
      out[[nm]] <- tt
      next
    }
    
    fdat <- Biobase::fData(eset)
    if (nrow(fdat) == 0L) {
      warning("  -> fData vacío para ", nm, ". No puedo mapear símbolos.")
      out[[nm]] <- tt
      next
    }
    
    # Find candidate symbol column
    cand <- symbol_candidates[symbol_candidates %in% colnames(fdat)]
    if (length(cand) == 0L) {
      warning("  -> No encontré ninguna columna de símbolo en fData de ", nm,
              ". Candidatas probadas: ",
              paste(symbol_candidates, collapse = ", "))
      out[[nm]] <- tt
      next
    }
    
    sym_col <- cand[1]
    message("  -> Usando columna de símbolo '", sym_col, "' de fData.")
    
    # Extract symbol map
    sym_map <- fdat[, sym_col, drop = TRUE]
    
    # Align with probes in the limma table
    common_probes <- intersect(rownames(tt), names(sym_map))
    if (length(common_probes) == 0L) {
      warning("  -> No hubo intersección entre rownames(tt) y rownames(fData) en ", nm)
      out[[nm]] <- tt
      next
    }
    
    # Create symbol vector aligned to the order of tt
    gene_sym <- rep(NA_character_, nrow(tt))
    names(gene_sym) <- rownames(tt)
    gene_sym[common_probes] <- sym_map[common_probes]
    
    tt[[out_col]] <- gene_sym
    
    out[[nm]] <- tt
  }
  
  return(out)
}

## ============================================================
## 6) Collapse by gene symbol using B statistic
## ============================================================

collapse_limma_table <- function(tt,
                                 symbol_col   = "Gene.symbol",
                                 keep_na_sym  = FALSE,
                                 split_multi  = FALSE,
                                 multi_sep    = " /// ") {
  stopifnot(is.data.frame(tt))
  stopifnot(symbol_col %in% colnames(tt))
  if (!"B" %in% colnames(tt)) {
    stop("La tabla no tiene columna 'B'; ¿seguro viene de limma::topTable con 'B'?")
  }
  
  df <- tt
  
  # Optional: split multiple symbols "A /// B /// C"
  if (split_multi) {
    # Requires tidyr/dplyr; can be reimplemented in base R if needed
    df <- tidyr::separate_rows(df,
                               !!rlang::sym(symbol_col),
                               sep = multi_sep)
  }
  
  # Normalise symbol as character
  df[[symbol_col]] <- as.character(df[[symbol_col]])
  
  if (!keep_na_sym) {
    df <- df[!is.na(df[[symbol_col]]) & df[[symbol_col]] != "", , drop = FALSE]
  }
  
  # If no rows remain after filtering, return an empty data.frame
  if (nrow(df) == 0L) {
    warning("collapse_limma_table: no quedaron filas con símbolo de gen válido.")
    return(df)
  }
  
  # Collapse by symbol, keeping the row with the highest B statistic
  # (ties are broken arbitrarily)
  o <- order(df[[symbol_col]], -df[["B"]])
  df_sorted <- df[o, , drop = FALSE]
  
  # Keep the first occurrence of each symbol (highest B)
  keep_idx <- !duplicated(df_sorted[[symbol_col]])
  df_collapsed <- df_sorted[keep_idx, , drop = FALSE]
  
  rownames(df_collapsed) <- df_collapsed[[symbol_col]]
  
  df_collapsed
}


collapse_limma_list <- function(de_list,
                                symbol_col   = "Gene.symbol",
                                keep_na_sym  = FALSE,
                                split_multi  = FALSE,
                                multi_sep    = " /// ") {
  stopifnot(is.list(de_list))
  
  out <- list()
  for (nm in names(de_list)) {
    message("Colapsando ", nm, " por ", symbol_col, " (max B)...")
    tt <- de_list[[nm]]
    if (is.null(tt)) {
      warning("  -> ", nm, " es NULL, lo omito.")
      next
    }
    out[[nm]] <- collapse_limma_table(
      tt,
      symbol_col  = symbol_col,
      keep_na_sym = keep_na_sym,
      split_multi = split_multi,
      multi_sep   = multi_sep
    )
  }
  out
}


## ============================================================
## 8) Robust Rank Aggregation
## ============================================================

make_rank_lists <- function(de_list,
                            direction = c("up", "down"),
                            p_col     = "adj.P.Val",
                            lfc_col   = "logFC",
                            stat_col  = NULL,     # si ya tienes un "stat" úsalo aquí
                            max_genes = Inf) {

  direction <- match.arg(direction)
  out <- list()

  for (nm in names(de_list)) {
    df <- de_list[[nm]]
    if (is.null(df) || nrow(df) == 0L) next

    # valid symbols
    if (!("Gene.symbol" %in% colnames(df))) {
      warning("En ", nm, " falta Gene.symbol. Omito.")
      next
    }
    df <- df[!is.na(df$Gene.symbol) & df$Gene.symbol != "", , drop = FALSE]
    if (nrow(df) == 0L) next

    # direction by logFC (used only to separate up/down)
    if (!lfc_col %in% colnames(df)) {
      warning("En ", nm, " falta ", lfc_col, ". Omito.")
      next
    }
    if (direction == "up") {
      df <- df[df[[lfc_col]] > 0, , drop = FALSE]
    } else {
      df <- df[df[[lfc_col]] < 0, , drop = FALSE]
    }
    if (nrow(df) == 0L) next

    # continuous ranking
    if (!is.null(stat_col) && stat_col %in% colnames(df)) {
      # higher value = stronger evidence
      df <- df[order(-df[[stat_col]]), , drop = FALSE]
    } else {
      # fallback: rank by p first, then |logFC|
      if (!p_col %in% colnames(df)) {
        warning("En ", nm, " falta ", p_col, " y no diste stat_col. Omito.")
        next
      }
      df <- df[order(df[[p_col]], -abs(df[[lfc_col]])), , drop = FALSE]
    }

    if (is.finite(max_genes) && nrow(df) > max_genes) {
      df <- df[seq_len(max_genes), , drop = FALSE]
    }

    out[[nm]] <- df$Gene.symbol
  }

  out <- out[sapply(out, length) > 0]
  out
}


## ============================================================
## Execution
## ============================================================
# 1) IDs de los datasets del paper de ovario
gse_ids <- c("GSE14407", "GSE18520", "GSE27651",
             "GSE38666", "GSE40595", "GSE54388")

esets <- download_gse_series(
  gse_ids,
  destdir  = "geo_data",
  platform = "GPL570",
  retries  = 2
)

# 2) Inspect phenoData to decide columns/labels:
cond_info <- inspect_gse_conditions(esets, max_levels = 12)

# 3) Define contrast_specs manually based on cond_info output:

contrast_specs <- list(
  GSE14407 = list(
    group_col = "disease state:ch1",
    case      = c("ovarian adenocarcinoma"),
    control   = c("normal")
  ),
  GSE18520 = list(
    group_col = "source_name_ch1",
    case      = c("papillary serous ovarian adenocarcinoma"),
    control   = c("normal ovarian surface epithelium (OSE)")
  ),
  GSE27651 = list(
    group_col = "description",
    case      = "leftover",
    control   = c("normal human ovarian surface epithelials cells")
  ),
  GSE38666 = list(
    group_col = "grade:ch1",
    case      = "leftover",
    control   = NA_character_
  ),
  GSE40595 = list(
    group_col = "source_name_ch1",
    case      = "leftover",
    control   = c(
    	"Microdissected normal ovarian stroma",
    	"Microdissected ovarian surface epthelium"
  	)
  ),
  GSE54388 = list(
    group_col = "diagnosis:ch1",
    case      = c("serous ovarian cancer"),
    control   = c("healthy (normal)")
  )
)

# 4) Check contrast_specs:
diagnostics <- check_contrast_specs(esets, contrast_specs)

# 5) Run limma on all datasets:
de_results <- limma_de_all(esets, contrast_specs)

# 6) Add Gene.symbol column using gconvert
de_results_sym <- add_gene_symbol_with_gconvert_safe(
  de_list = de_results,
  id_col  = "ID",
  organism = "hsapiens"
)


# 7) Collapse by symbol using the function defined above
de_genes <- collapse_limma_list(
  de_list     = de_results_sym,
  symbol_col  = "Gene.symbol",
  keep_na_sym = FALSE,
  split_multi = FALSE
)
# Example: inspect a dataset's header
head(de_genes$GSE27651)
saveRDS(de_genes, "DEGs_alldsets_GEO.rds")

# 8) Sort lists for RRA
rank_lists_up <- make_rank_lists(de_genes, direction="up",   p_col="adj.P.Val")
rank_lists_dn <- make_rank_lists(de_genes, direction="down", p_col="adj.P.Val")

# 9) Run RRA
# UP-regulated meta-set
rra_up <- RobustRankAggreg::aggregateRanks(
  glist = rank_lists_up,
  full  = TRUE   # para tener ranking en todo el universo posible
)

# DOWN-regulated meta-set
rra_down <- RobustRankAggreg::aggregateRanks(
  glist = rank_lists_dn,
  full  = TRUE
)

adjust_rra_no_filter <- function(rra_tbl, method = "BH") {
  stopifnot(all(c("Name", "Score") %in% colnames(rra_tbl)))
  
  out <- data.frame(
    Gene.symbol = as.character(rra_tbl$Name),
    p_raw       = as.numeric(rra_tbl$Score),
    stringsAsFactors = FALSE
  )
  
  out$p_adj <- p.adjust(out$p_raw, method = method)
  out <- out[order(out$p_adj, out$p_raw), , drop = FALSE]
  out$rank_adj <- rank(out$p_adj, ties.method = "average")
  out
}

# merge_rra_up_down <- function(rra_up, rra_down, method = "BH") {
#   up   <- adjust_rra_no_filter(rra_up, method = method)
#   down <- adjust_rra_no_filter(rra_down, method = method)
  
#   up$direction   <- "up"
#   down$direction <- "down"
  
#   geo_rra <- rbind(up, down)
#   geo_rra <- geo_rra[order(geo_rra$rank_adj, geo_rra$p_adj, geo_rra$p_raw), ]
  
#   rownames(geo_rra) <- NULL
#   geo_rra
# }

# geo_rra_all <- merge_rra_up_down(rra_up, rra_down, method = "BH")

library(dplyr)

# =========================
# 1) Adjusted UP/DOWN
# =========================
up <- adjust_rra_no_filter(rra_up, method = "BH") %>%
  transmute(
    gene     = as.character(Gene.symbol),
    p_raw_up = as.numeric(p_raw),
    p_adj_up = as.numeric(p_adj)
  )

down <- adjust_rra_no_filter(rra_down, method = "BH") %>%
  transmute(
    gene     = as.character(Gene.symbol),
    p_raw_dn = as.numeric(p_raw),
    p_adj_dn = as.numeric(p_adj)
  )

# =========================
# 2) m = full universe (without collapsing)
# =========================
m <- full_join(up, down, by = "gene") %>%
  filter(!is.na(gene), nzchar(gene))

# =========================
# 3) geo_rra_1row = 1 row per gene (direction chosen by smallest p_adj)
# =========================
geo_rra_1row <- m %>%
  mutate(
    direction = case_when(
      is.na(p_adj_dn) ~ "up",
      is.na(p_adj_up) ~ "down",
      p_adj_up <= p_adj_dn ~ "up",
      TRUE ~ "down"
    ),
    p_adj = pmin(p_adj_up, p_adj_dn, na.rm = TRUE),
    p_raw = pmin(p_raw_up, p_raw_dn, na.rm = TRUE)
  ) %>%
  select(gene, direction, p_raw, p_adj) %>%
  filter(is.finite(p_adj)) %>%
  arrange(p_adj, p_raw) %>%
  mutate(
    rank_adj = rank(p_adj, ties.method = "average"),
    rank_raw = rank(p_raw, ties.method = "average")
  ) %>%
  group_by(direction) %>%
  arrange(p_adj, p_raw, .by_group = TRUE) %>%
  mutate(rank_adj_within_dir = row_number()) %>%
  ungroup()


head(geo_rra_1row)
table(geo_rra_1row$direction)


saveRDS(geo_rra_1row, file = "GEO_ovarian_cancer_RRA_1row.rds")

# =========================
# 4) Tabla "reportable" genome-wide sin filtros
#    (rank para todos los genes, con score Z firmado)
# =========================

rra_report_all <- geo_rra_1row %>%
  mutate(
    sign_dir = ifelse(direction == "up", 1, -1),

    # clamp01 inline: evitar 0 y 1 exactos
    p_adj_clamped = pmin(
      pmax(as.numeric(p_adj), 1e-300),
      1 - 1e-300
    ),

    z_abs = qnorm(1 - p_adj_clamped / 2),
    rra_z = sign_dir * z_abs,
    abs_rra_z = abs(rra_z)
  ) %>%
  arrange(desc(abs_rra_z), p_adj, p_raw) %>%
  mutate(rank_by_absZ = row_number()) %>%
  group_by(direction) %>%
  mutate(rank_by_absZ_within_dir = row_number()) %>%
  ungroup() %>%
  select(
    gene, direction,
    p_raw, p_adj,
    rank_adj, rank_raw, rank_adj_within_dir,
    rra_z, abs_rra_z,
    rank_by_absZ, rank_by_absZ_within_dir
  )

# Column descriptions for rra_report_all
# gene: HGNC gene symbol. One row per gene, genome-wide.
# direction: RRA effect direction (“up” or “down”), chosen as the side (UP/DOWN) with the smallest p_adj.
# p_raw: Raw RRA p-value for the winning direction (UP or DOWN).
# p_adj: BH-adjusted RRA p-value for the winning direction. This is the formal p-value reported in a standard RRA table.

# rank_adj: Global rank (all directions) by ascending p_adj. Ties resolved with ties.method = “average”.
# rank_raw: Global rank by ascending p_raw. Useful only for inspection.
# rank_adj_within_dir: Rank by p_adj within each direction (up and down separately). Equivalent to “Top up genes” and “Top down genes” in classic tables.

# rra_z: Signed continuous RRA score. Magnitude = evidence strength; sign = direction.
# abs_rra_z: Absolute value of rra_z. Measures RRA signal strength regardless of direction.

# * rank_by_absZ: Final global rank by RRA evidence, ordered by descending abs_rra_z. This is the primary genome-wide ranking.
# * rank_by_absZ_within_dir: Rank by RRA evidence within each direction (up / down). Useful for reporting “Top upregulated” and “Top downregulated” genes.


saveRDS(rra_report_all, file = "GEO_ovarian_cancer_RRA_results.rds")
write.table(rra_report_all,
								file = "GEO_ovarian_cancer_RRA_results.tsv",
								sep  = "\t",
								row.names = FALSE,
								quote = FALSE)



# =========================
# CONTROL: Random subsampling to rule out sample-size bias
# =========================

# Verifies that the correlation improvement when restricting to genes with
# strong RRA signal is not a sample-size artefact but a genuine property
# of the most reproducible genes.

library(readr)

# --- Load DESeq2 statistics from TCGA ---
DE_AE   <- readRDS("../1_DGE_AE/DE_full_OVARY_DESeq2_AE.rds")
DE_GTEx <- readRDS("../0_DGE_GTEx/DE_full_OVARY_DESeq2_GTEx.rds")

# --- Deduplicate per gene (abs_max) before joining ---
de_ae_uniq <- DE_AE %>%
  select(gene = Symbol, stat_ae = stat) %>%
  mutate(gene = as.character(gene)) %>%
  filter(!is.na(gene), nzchar(gene), is.finite(stat_ae)) %>%
  group_by(gene) %>%
  slice_max(abs(stat_ae), n = 1, with_ties = FALSE) %>%
  ungroup()

de_gtex_uniq <- DE_GTEx %>%
  select(gene = Symbol, stat_gtex = stat) %>%
  mutate(gene = as.character(gene)) %>%
  filter(!is.na(gene), nzchar(gene), is.finite(stat_gtex)) %>%
  group_by(gene) %>%
  slice_max(abs(stat_gtex), n = 1, with_ties = FALSE) %>%
  ungroup()

# --- Build joint table ---
df_full <- rra_report_all %>%
  inner_join(de_ae_uniq,   by = "gene") %>%
  inner_join(de_gtex_uniq, by = "gene") %>%
  filter(is.finite(rra_z), is.finite(stat_ae), is.finite(stat_gtex))

n_total <- nrow(df_full)

# --- Genome-wide correlation (baseline) ---
rho_ae_full   <- cor(df_full$rra_z, df_full$stat_ae,   method = "spearman")
rho_gtex_full <- cor(df_full$rra_z, df_full$stat_gtex, method = "spearman")
message("Correlación genome-wide — AE: ", round(rho_ae_full, 3),
        " | GTEx: ", round(rho_gtex_full, 3))

# --- Correlation in the subset of genes with strong RRA signal ---
# Uses the same criterion as in the manuscript (top N by abs_rra_z)
n_sub <- 2000  # ajusta al valor real que usaste
df_top <- df_full %>% arrange(desc(abs_rra_z)) %>% slice_head(n = n_sub)

rho_ae_top   <- cor(df_top$rra_z, df_top$stat_ae,   method = "spearman")
rho_gtex_top <- cor(df_top$rra_z, df_top$stat_gtex, method = "spearman")
message("Correlación top-", n_sub, " RRA — AE: ", round(rho_ae_top, 3),
        " | GTEx: ", round(rho_gtex_top, 3))

# --- Null distribution by random subsampling (robustness check) ---
set.seed(42)
n_perm <- 1000

null_cors <- replicate(n_perm, {
  idx <- sample(n_total, n_sub)
  df_s <- df_full[idx, ]
  c(
    ae   = cor(df_s$rra_z, df_s$stat_ae,   method = "spearman"),
    gtex = cor(df_s$rra_z, df_s$stat_gtex, method = "spearman")
  )
})

null_ae   <- null_cors["ae",   ]
null_gtex <- null_cors["gtex", ]

# empirical p-value (proportion of nulls >= observed)
p_ae   <- mean(abs(null_ae)   >= abs(rho_ae_top))
p_gtex <- mean(abs(null_gtex) >= abs(rho_gtex_top))

message("p-valor empírico (submuestreo) — AE: ", p_ae, " | GTEx: ", p_gtex)

# --- Supplementary figure ---
library(ggplot2)

df_null <- data.frame(
  rho    = c(null_ae, null_gtex),
  source = rep(c("AE", "GTEx"), each = n_perm)
)
obs_df <- data.frame(
  rho    = c(rho_ae_top, rho_gtex_top),
  source = c("AE", "GTEx")
)

p_joel <- ggplot(df_null, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(data = obs_df, aes(xintercept = rho),
             color = "red3", linewidth = 1, linetype = "dashed") +
  facet_wrap(~source, ncol = 2) +
  labs(
    title = paste0("Null distribution: random subsets of n = ", n_sub, " genes"),
    subtitle = "Red dashed line = observed correlation in top-RRA subset",
    x = "Spearman ρ (RRA z-score vs DESeq2 stat)",
    y = "Count"
  ) +
  theme_bw(base_size = 13)

ggsave("FigS_submuestreo_control_Joel.png", p_joel, width = 9, height = 4, dpi = 300)
ggsave("FigS_submuestreo_control_Joel.pdf", p_joel, width = 9, height = 4)