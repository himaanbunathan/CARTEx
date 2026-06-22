# Internal helpers for single-cell CARTEx scoring.
# None of these are exported except the S3 print method for the report object.

#' @keywords internal
#' @noRd
.cartex_gene_map <- function() {
  path <- system.file("extdata", "cartex_gene_map.rds", package = "CARTEx")
  if (path == "") {
    # Fall back to the bare weights vector if the richer map is unavailable.
    wpath <- system.file("extdata", "weights.rds", package = "CARTEx")
    if (wpath == "") stop("CARTEx gene weights not found in the installed package.")
    w <- readRDS(wpath)
    return(data.frame(
      cartex_symbol = names(w),
      clean_symbol  = gsub("\\.", "-", names(w)),
      ensembl_id    = NA_character_,
      weight        = as.numeric(w),
      stringsAsFactors = FALSE
    ))
  }
  readRDS(path)
}

#' @keywords internal
#' @noRd
.zscore <- function(s) {
  s <- as.numeric(s)
  (s - mean(s, na.rm = TRUE)) / stats::sd(s, na.rm = TRUE)
}

#' Match CARTEx genes against an object's features.
#'
#' Tries, in order: exact symbol, cleaned symbol ("." -> "-" / curated renames),
#' Ensembl ID against supplied feature IDs, then Ensembl ID against the feature
#' names themselves (objects that use Ensembl IDs as row names). Records which
#' rule matched each gene so the decision is fully auditable.
#'
#' @param features Character vector of feature names (gene symbols / row names).
#' @param feature_ids Optional character vector of Ensembl IDs, parallel to
#'   `features` (e.g. AnnData `var/_index`).
#' @param gene_map The CARTEx gene map data.frame.
#' @return A list: `matched` (data.frame), `unmatched` (character), `all`,
#'   `n_requested`, `n_matched`.
#' @keywords internal
#' @noRd
.match_cartex_genes <- function(features, feature_ids = NULL,
                                gene_map = .cartex_gene_map()) {
  features <- as.character(features)
  keep_f <- !duplicated(features)
  sym_lookup <- stats::setNames(seq_along(features)[keep_f], features[keep_f])

  id_lookup <- NULL
  if (!is.null(feature_ids)) {
    feature_ids <- as.character(feature_ids)
    keep_i <- !duplicated(feature_ids)
    id_lookup <- stats::setNames(seq_along(feature_ids)[keep_i], feature_ids[keep_i])
  }

  look <- function(tbl, key) {
    if (is.null(tbl) || is.na(key)) return(NA_integer_)
    v <- tbl[key]
    if (is.na(v)) NA_integer_ else as.integer(v)
  }

  n <- nrow(gene_map)
  idx <- integer(n); rule <- character(n); matched <- character(n)
  for (k in seq_len(n)) {
    cs <- gene_map$cartex_symbol[k]
    ks <- gene_map$clean_symbol[k]
    es <- gene_map$ensembl_id[k]
    i <- look(sym_lookup, cs); r <- "symbol_exact"; mv <- cs
    if (is.na(i) && !identical(ks, cs)) { i <- look(sym_lookup, ks); r <- "symbol_clean"; mv <- ks }
    if (is.na(i)) { i <- look(id_lookup, es); r <- "ensembl_id"; mv <- es }
    if (is.na(i)) { i <- look(sym_lookup, es); r <- "ensembl_in_features"; mv <- es }
    idx[k] <- i
    rule[k] <- if (is.na(i)) NA_character_ else r
    matched[k] <- if (is.na(i)) NA_character_ else mv
  }

  all <- data.frame(
    cartex_symbol = gene_map$cartex_symbol,
    weight = gene_map$weight,
    feature_index = idx,
    match_rule = rule,
    matched_feature = matched,
    stringsAsFactors = FALSE
  )
  ok <- !is.na(all$feature_index)
  list(
    matched = all[ok, , drop = FALSE],
    unmatched = all$cartex_symbol[!ok],
    all = all,
    n_requested = n,
    n_matched = sum(ok)
  )
}

#' Classify an expression matrix as raw counts vs (log-)normalized.
#'
#' Value-based heuristics that deliberately override any declared label, because
#' AnnData `uns/X_name` is frequently inaccurate.
#'
#' @param values Numeric vector of (preferably nonzero) sampled expression values.
#' @param declared Optional declared layer name/label (used only in the reason).
#' @keywords internal
#' @noRd
.inspect_expression <- function(values, declared = NULL) {
  v <- values[is.finite(values)]
  v <- v[v != 0]
  if (!length(v)) {
    return(list(state = "unknown", is_integer = NA, max_value = NA_real_,
                reason = "no nonzero values sampled"))
  }
  is_int <- all(abs(v - round(v)) < 1e-6)
  mx <- max(v)
  state <- if (is_int) {
    "counts"
  } else if (mx <= 30) {
    "log_normalized"
  } else {
    "linear_normalized"
  }
  reason <- sprintf(
    "%s values; max=%.3g; %s%s",
    if (is_int) "integer" else "non-integer", mx,
    switch(state,
      counts = "integers => raw counts",
      log_normalized = "non-integer, small range => log-normalized",
      linear_normalized = "non-integer, large range => linearly normalized (not logged)"),
    if (!is.null(declared)) sprintf(" [declared: '%s']", declared) else ""
  )
  list(state = state, is_integer = is_int, max_value = mx, reason = reason)
}

#' Apply normalization to bring a matrix to log-normalized space.
#'
#' @param expr genes x cells matrix (dense or Matrix) of the *selected* layer.
#' @param state One of "counts", "linear_normalized", "log_normalized".
#' @param cell_total Per-cell total library size (sum over ALL genes) of the
#'   selected layer; required when `state == "counts"`.
#' @param scale_factor LogNormalize scale factor (default 1e4).
#' @return list(expr = normalized matrix, action = character).
#' @keywords internal
#' @noRd
.normalize_expr <- function(expr, state, cell_total = NULL, scale_factor = 1e4) {
  if (state == "log_normalized") {
    return(list(expr = expr, action = "none (already log-normalized)"))
  }
  if (state == "linear_normalized") {
    return(list(expr = log1p(expr),
                action = "log1p (linearly-normalized input)"))
  }
  if (state == "counts") {
    if (is.null(cell_total)) {
      stop("Library sizes (cell_total) are required to LogNormalize raw counts.")
    }
    cell_total[cell_total == 0] <- 1
    # value' = log1p(count / total * scale_factor), applied per cell (column).
    cpm <- Matrix::t(Matrix::t(expr) / cell_total) * scale_factor
    return(list(expr = log1p(cpm),
                action = sprintf("LogNormalize(scale.factor=%g)", scale_factor)))
  }
  list(expr = expr, action = "none (unknown state; left unchanged)")
}

#' Per-cell CARTEx score from a matched-gene expression matrix.
#'
#' @param expr genes x cells matrix; row order must align to `weights`.
#' @param weights Numeric weight per row (gene).
#' @keywords internal
#' @noRd
.score_cells <- function(expr, weights) {
  if (nrow(expr) != length(weights)) {
    stop("Expression rows and weights length differ.")
  }
  raw <- as.numeric(Matrix::crossprod(expr, matrix(as.numeric(weights), ncol = 1)))
  names(raw) <- colnames(expr)
  z <- .zscore(raw)
  data.frame(
    cartex_score = z,
    cartex_score_raw = raw,
    cartex_score_int = .integerize(z),
    row.names = colnames(expr),
    stringsAsFactors = FALSE
  )
}

#' Detect common QC metadata columns and summarize them.
#'
#' @param meta Cell metadata data.frame.
#' @keywords internal
#' @noRd
.detect_qc <- function(meta) {
  nm <- names(meta)
  find <- function(pattern) {
    hits <- nm[grepl(pattern, nm, ignore.case = TRUE)]
    hits[vapply(hits, function(h) is.numeric(meta[[h]]), logical(1))]
  }
  patterns <- list(
    n_features = "^(nFeature_RNA|n_genes_by_counts|n_genes|nGene|detected_genes)$",
    n_counts   = "^(nCount_RNA|n_counts|total_counts|nUMI|library_size)$",
    pct_mito   = "(percent\\.?mt|pct_counts_mt|percent_?mito|mito.*(percent|frac)|^mt_frac)"
  )
  detected <- list()
  for (key in names(patterns)) {
    cols <- find(patterns[[key]])
    if (length(cols)) {
      col <- cols[1]
      x <- meta[[col]]
      detected[[key]] <- list(
        column = col,
        min = min(x, na.rm = TRUE),
        median = stats::median(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE)
      )
    }
  }
  # Heuristic: does it look like QC filtering was already applied?
  filtered <- NA
  reasons <- character(0)
  if (!is.null(detected$n_features)) {
    if (detected$n_features$min >= 100) {
      filtered <- TRUE
      reasons <- c(reasons, sprintf("min %s = %g (>=100)",
                                    detected$n_features$column, detected$n_features$min))
    } else {
      filtered <- FALSE
      reasons <- c(reasons, sprintf("min %s = %g (<100; cells with few genes remain)",
                                    detected$n_features$column, detected$n_features$min))
    }
  }
  if (!is.null(detected$pct_mito)) {
    if (detected$pct_mito$max <= 25) {
      reasons <- c(reasons, sprintf("max %s = %.1f (<=25%%)",
                                    detected$pct_mito$column, detected$pct_mito$max))
      if (is.na(filtered)) filtered <- TRUE
    } else {
      reasons <- c(reasons, sprintf("max %s = %.1f (>25%%; high-mito cells remain)",
                                    detected$pct_mito$column, detected$pct_mito$max))
    }
  }
  list(detected = detected,
       filtering_likely_applied = filtered,
       filtering_reasons = reasons)
}
