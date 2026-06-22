#' Compute the CARTEx exhaustion score for single-cell data
#'
#' Scores single-cell data with the published CARTEx CD8+ T-cell exhaustion
#' signature (200 PCA-weighted genes). Built to run on a laptop even for large
#' files: for HDF5-backed formats it reads only the CARTEx genes (plus a few CD8
#' marker genes) from disk rather than loading the whole matrix.
#'
#' Because CARTEx is only valid for CD8+ T cells, every cell is classified as
#' CD8+ or not before scoring (see `cd8`/`cd8_method`). And because a per-cell
#' (or per-cluster) score inherits choices baked into the input object, every
#' decision the function makes is captured in a `cartex_report` object so the
#' assumptions behind the score are explicit and auditable.
#'
#' @section Analysis level:
#' With `level = "single_cell"` (default) each cell is scored individually. With
#' `level = "pseudobulk"`, cells are grouped (by `group_by`, e.g. cluster), their
#' raw counts are summed and LogNormalized, and the resulting pseudobulk profiles
#' are scored -- following the reference CARTEx pipeline
#' (Seurat::AggregateExpression-style aggregation). Set `pseudo_replicates > 1`
#' to randomly split each group into replicate pseudobulk samples.
#'
#' @section Score definition:
#' For each cell or pseudobulk sample, the raw score is the dot product of its
#' log-normalized expression with the 200 CARTEx weights (`cartex_score_raw`).
#' `cartex_score` is that value Z-scored across all scored units, and
#' `cartex_score_int` is the Z-score rounded and capped at +/-4 (matching the
#' reference `integerize`). The Z-score is population-relative, so scores are
#' comparable only within one scored object.
#'
#' @param input A file path (`.rds`, `.RData`, `.h5ad`, `.h5Seurat`, `.loom`),
#'   an in-memory Seurat object, or a gene x cell matrix.
#' @param level `"single_cell"` (default) or `"pseudobulk"`.
#' @param cd8 How to use the CD8+ classification: `"flag"` (default) scores all
#'   cells and records a `cd8_status` column; `"restrict"` scores only confident
#'   CD8+ T cells.
#' @param cd8_method `"auto"` (default: cell-type annotation if present, else
#'   marker genes), `"metadata"`, `"markers"`, or `"none"`.
#' @param cd8_column Optional explicit cell-type annotation column name.
#' @param group_by For pseudobulk: metadata column(s) to aggregate by. If `NULL`,
#'   a cluster-like column is auto-detected.
#' @param pseudo_replicates For pseudobulk: number of random replicate bins per
#'   group (default 1 = none).
#' @param seed Random seed for pseudo-replicate assignment (default 1).
#' @param normalization `"auto"` (detect counts vs normalized), `"lognorm"`, or
#'   `"assume_normalized"`.
#' @param assay,layer Optional assay / layer (Seurat / h5Seurat).
#' @param min_match_frac Warn if fewer than this fraction of CARTEx genes match.
#' @param return `"object"` (default), `"scores"`, or `"report"`.
#' @param verbose Print the report summary when done (default TRUE).
#' @return Depending on `return` and `level`: the annotated object, a scores
#'   data.frame, or the `cartex_report`.
#' @seealso [cartex_score()] for bulk data, [read_single_cell()] for loading only.
#' @export
cartex_sc <- function(input,
                      level = c("single_cell", "pseudobulk"),
                      cd8 = c("flag", "restrict"),
                      cd8_method = c("auto", "metadata", "markers", "none"),
                      cd8_column = NULL,
                      group_by = NULL,
                      pseudo_replicates = 1L,
                      seed = 1L,
                      normalization = c("auto", "lognorm", "assume_normalized"),
                      assay = NULL, layer = NULL,
                      min_match_frac = 0.5,
                      return = c("object", "scores", "report"),
                      verbose = TRUE) {
  level <- match.arg(level)
  cd8 <- match.arg(cd8)
  cd8_method <- match.arg(cd8_method)
  normalization <- match.arg(normalization)
  return <- match.arg(return)
  t0 <- Sys.time()
  warnings_acc <- character(0)

  gene_map <- .cartex_gene_map()
  want <- if (level == "pseudobulk") "counts" else "normalized"
  loaded <- .load_expression(input, gene_map, want = want, assay = assay, layer = layer)

  mm <- loaded$matched
  match_frac <- mm$n_matched / mm$n_requested
  if (mm$n_matched == 0)
    stop("No CARTEx genes matched the object's features. Check that genes are ",
         "labelled with HGNC symbols or Ensembl IDs.")
  if (match_frac < min_match_frac)
    warnings_acc <- c(warnings_acc, sprintf(
      "Only %d/%d CARTEx genes matched (%.0f%%); score may be unreliable.",
      mm$n_matched, mm$n_requested, 100 * match_frac))

  cartex_expr <- loaded$cartex_expr
  meta <- loaded$meta
  cells <- colnames(cartex_expr)
  # Align metadata rows to the expression columns.
  if (nrow(meta) && all(cells %in% rownames(meta))) meta <- meta[cells, , drop = FALSE]

  # --- CD8+ T-cell classification -------------------------------------------
  cd8res <- .classify_cd8(meta, marker_expr = loaded$marker_expr,
                          method = cd8_method, cd8_column = cd8_column)
  cd8_status <- cd8res$status[cells]

  keep <- rep(TRUE, length(cells))
  if (cd8 == "restrict") {
    keep <- as.character(cd8_status) == "CD8_pos"
    if (sum(keep) == 0)
      stop("cd8 = 'restrict' but no cells were classified as CD8+ T cells ",
           "(method: ", cd8res$method_used, "). Check the annotation/markers.")
    cartex_expr <- cartex_expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
    if (!is.null(loaded$cell_total)) loaded$cell_total <- loaded$cell_total[keep]
    cells <- cells[keep]
    warnings_acc <- c(warnings_acc, sprintf(
      "Restricted to %d CD8+ T cells of %d (%s).",
      sum(keep), length(keep), cd8res$method_used))
  }
  cell_total <- loaded$cell_total

  # --- Score: single-cell or pseudobulk -------------------------------------
  weights <- mm$matched$weight
  detected_state <- loaded$state %||% "unknown"

  if (level == "single_cell") {
    norm <- .apply_norm(cartex_expr, detected_state, normalization, cell_total)
    scores <- .score_cells(norm$expr, weights)
    scores$cd8_status <- as.character(cd8_status[rownames(scores)])
    zero_cells <- sum(Matrix::colSums(abs(norm$expr)) == 0)
    if (zero_cells > 0)
      warnings_acc <- c(warnings_acc, sprintf(
        "%d cell(s) have zero expression across all matched CARTEx genes.",
        zero_cells))
    norm_action <- norm$action
    pb_info <- NULL
  } else {
    grp <- .resolve_group_by(meta, group_by)
    keys <- .build_group_keys(meta, grp$columns, pseudo_replicates, seed)
    is_counts <- detected_state == "counts" && normalization != "assume_normalized"
    if (!is_counts)
      warnings_acc <- c(warnings_acc,
        "Pseudobulk on non-count data: averaging normalized expression instead ",
        "of summing counts.")
    pb <- .score_pseudobulk(cartex_expr, cell_total, keys, weights,
                            is_counts = is_counts)
    scores <- pb$scores
    norm_action <- pb$action
    zero_cells <- 0
    pb_info <- list(columns = grp$columns, source = grp$source,
                    n_pseudobulk = pb$n_pseudobulk,
                    replicates = pseudo_replicates, group = pb$group)
  }

  # --- Report ----------------------------------------------------------------
  qc <- .detect_qc(meta)
  report <- .build_report(
    provenance = loaded$provenance,
    level = level,
    dims = list(n_cells = loaded$n_cells, n_genes_total = loaded$n_genes,
                n_scored = nrow(scores)),
    matching = mm,
    normalization = list(detected = detected_state, requested = normalization,
                         action = norm_action, layer = loaded$layer_used,
                         declared = loaded$declared, reason = loaded$inspect_reason),
    cd8 = list(mode = cd8, method = cd8res$method_used, column = cd8res$column,
               counts = cd8res$counts, notes = cd8res$notes, kept = sum(keep)),
    pseudobulk = pb_info,
    qc = qc, meta_columns = names(meta),
    scores = scores, zero_cells = zero_cells, warnings = warnings_acc,
    elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  if (verbose) print(report)
  for (w in unique(warnings_acc)) warning(w, call. = FALSE)

  if (return == "report") return(report)
  if (return == "scores") { attr(scores, "cartex_report") <- report; return(scores) }

  .assemble_output(loaded, scores, cd8_status, level, pb_info, meta, cells, report)
}

#' Decide and apply normalization for single-cell scoring.
#' @keywords internal
#' @noRd
.apply_norm <- function(expr, detected_state, normalization, cell_total) {
  if (normalization == "assume_normalized")
    return(list(expr = expr, action = "none (assume_normalized requested)"))
  if (normalization == "lognorm") {
    if (is.null(cell_total))
      return(list(expr = log1p(expr), action = "log1p (no library sizes)"))
    return(.normalize_expr(expr, "counts", cell_total))
  }
  .normalize_expr(expr, detected_state, cell_total)
}

#' Auto-detect or validate a pseudobulk grouping column.
#' @keywords internal
#' @noRd
.resolve_group_by <- function(meta, group_by) {
  if (!is.null(group_by)) {
    miss <- setdiff(group_by, names(meta))
    if (length(miss))
      stop("group_by column(s) not found in metadata: ", paste(miss, collapse = ", "))
    return(list(columns = group_by, source = "user"))
  }
  candidates <- c("seurat_clusters", "clusters", "cluster", "leiden", "louvain",
                  "SCT_snn_res.0.8", "RNA_snn_res.0.8", "cell_type", "celltype",
                  "CellType", "niche", "annotation")
  hit <- intersect(candidates, names(meta))
  if (!length(hit))
    stop("No cluster-like column found for pseudobulk; specify group_by = ",
         "<metadata column name>. Available columns: ",
         paste(utils::head(names(meta), 20), collapse = ", "))
  list(columns = hit[1], source = "auto-detected")
}

#' Build per-cell pseudobulk group keys from grouping columns + replicate bins.
#' @keywords internal
#' @noRd
.build_group_keys <- function(meta, columns, pseudo_replicates, seed) {
  base <- do.call(paste, c(lapply(columns, function(c) as.character(meta[[c]])),
                           sep = " | "))
  if (pseudo_replicates > 1) {
    bins <- .pseudobulk_replicates(nrow(meta), pseudo_replicates, seed)
    base <- paste0(base, " | rep", bins)
  }
  base
}

#' Build the returned object (annotate Seurat / merge scores + metadata).
#' @keywords internal
#' @noRd
.assemble_output <- function(loaded, scores, cd8_status, level, pb_info, meta,
                             cells, report) {
  if (level == "single_cell") {
    if (!is.null(loaded$seurat)) {
      obj <- loaded$seurat
      obj$cartex_score     <- scores[colnames(obj), "cartex_score"]
      obj$cartex_score_raw <- scores[colnames(obj), "cartex_score_raw"]
      obj$cartex_score_int <- scores[colnames(obj), "cartex_score_int"]
      obj$cd8_status       <- as.character(cd8_status[colnames(obj)])
      attr(obj, "cartex_report") <- report
      return(obj)
    }
    out <- scores  # already includes cd8_status
    if (ncol(meta)) {
      common <- intersect(rownames(out), rownames(meta))
      out <- cbind(out[common, , drop = FALSE], meta[common, , drop = FALSE])
    }
    if (!is.null(loaded$umap)) {                 # carry UMAP coords for plotting
      u <- loaded$umap[rownames(out), , drop = FALSE]
      out$UMAP_1 <- u$UMAP_1; out$UMAP_2 <- u$UMAP_2
    }
    attr(out, "cartex_report") <- report
    return(out)
  }

  # pseudobulk: scores already carry pseudobulk_group + n_cells.
  out <- scores
  # When a live Seurat object exists, also map each group's score onto its cells
  # (handy for UMAP overlays), keeping the per-group table as an attribute.
  if (!is.null(loaded$seurat)) {
    obj <- loaded$seurat
    cell_keys <- .build_group_keys(obj@meta.data, pb_info$columns, 1L, 1L)
    obj$cartex_pseudobulk     <- scores[cell_keys, "cartex_score"]
    obj$cartex_pseudobulk_raw <- scores[cell_keys, "cartex_score_raw"]
    attr(obj, "cartex_report") <- report
    attr(obj, "cartex_pseudobulk_scores") <- out
    return(obj)
  }
  attr(out, "cartex_report") <- report
  out
}

#' Assemble the CARTEx transparency report.
#' @keywords internal
#' @noRd
.build_report <- function(provenance, level, dims, matching, normalization, cd8,
                          pseudobulk, qc, meta_columns, scores, zero_cells,
                          warnings, elapsed) {
  structure(list(
    provenance = provenance, level = level, dims = dims,
    gene_matching = list(
      n_requested = matching$n_requested, n_matched = matching$n_matched,
      match_fraction = matching$n_matched / matching$n_requested,
      by_rule = table(matching$matched$match_rule), unmatched = matching$unmatched),
    normalization = normalization, cd8 = cd8, pseudobulk = pseudobulk, qc = qc,
    metadata_columns = meta_columns,
    score_summary = summary(scores$cartex_score_raw), zero_cells = zero_cells,
    caveat = paste("cartex_score is Z-scored across the units in THIS object;",
                   "it is relative to this population and not comparable across",
                   "separately scored datasets."),
    warnings = unique(warnings), elapsed_secs = elapsed),
    class = "cartex_report")
}

#' Print a CARTEx transparency report
#' @param x A `cartex_report` object.
#' @param ... Ignored.
#' @export
print.cartex_report <- function(x, ...) {
  hr <- function() cat(strrep("-", 64), "\n", sep = "")
  cat("== CARTEx single-cell report ==\n"); hr()
  cat(sprintf("Input    : %s\n", x$provenance$format))
  if (is.character(x$provenance$source))
    cat(sprintf("Source   : %s\n", x$provenance$source))
  cat(sprintf("Level    : %s\n", x$level))
  cat(sprintf("Cells    : %s in object\n", format(x$dims$n_cells, big.mark = ",")))
  cat(sprintf("Scored   : %s %s\n", format(x$dims$n_scored, big.mark = ","),
              if (x$level == "pseudobulk") "pseudobulk samples" else "cells"))
  cat(sprintf("Layer    : %s\n", x$normalization$layer)); hr()

  gm <- x$gene_matching
  cat(sprintf("Gene match : %d / %d CARTEx genes (%.0f%%)\n",
              gm$n_matched, gm$n_requested, 100 * gm$match_fraction))
  if (length(gm$by_rule))
    cat("  by rule  :", paste(sprintf("%s=%d", names(gm$by_rule), gm$by_rule),
                              collapse = ", "), "\n")
  if (length(gm$unmatched))
    cat("  unmatched:", paste(gm$unmatched, collapse = ", "), "\n")
  hr()

  cd <- x$cd8
  cat(sprintf("CD8+ T cells: method '%s'%s, mode '%s'\n", cd$method,
              if (!is.na(cd$column)) sprintf(" [%s]", cd$column) else "", cd$mode))
  comp <- cd$counts[cd$counts > 0]
  if (length(comp))
    cat("  composition:", paste(sprintf("%s=%d", names(comp), comp),
                                collapse = ", "), "\n")
  for (n in cd$notes) cat("  note:", n, "\n")
  hr()

  nz <- x$normalization
  cat(sprintf("Normalization: detected '%s' -> %s\n", nz$detected, nz$action))
  if (!is.null(nz$reason)) cat(sprintf("  basis    : %s\n", nz$reason))
  hr()

  if (!is.null(x$pseudobulk)) {
    pb <- x$pseudobulk
    cat(sprintf("Pseudobulk : grouped by %s (%s) -> %d samples%s\n",
                paste(pb$columns, collapse = " + "), pb$source, pb$n_pseudobulk,
                if (pb$replicates > 1) sprintf(", %d replicates", pb$replicates) else ""))
    hr()
  }

  if (length(x$qc$detected)) {
    cat("QC fields detected:\n")
    for (k in names(x$qc$detected)) {
      d <- x$qc$detected[[k]]
      cat(sprintf("  %-10s [%s]: min=%.3g median=%.3g max=%.3g\n",
                  k, d$column, d$min, d$median, d$max))
    }
    if (!is.na(x$qc$filtering_likely_applied))
      cat(sprintf("  QC filtering appears %s\n",
                  if (x$qc$filtering_likely_applied) "ALREADY applied" else
                    "NOT applied (or permissive)"))
  } else {
    cat("QC fields detected: none (no nFeature/nCount/percent.mt columns found)\n")
  }
  cat(sprintf("Metadata columns: %d available\n", length(x$metadata_columns))); hr()

  cat("Raw score summary:\n"); print(x$score_summary)
  if (x$zero_cells > 0)
    cat(sprintf("Note: %d cell(s) had zero matched-gene signal.\n", x$zero_cells))
  cat("\nCaveat:", x$caveat, "\n")
  if (length(x$warnings)) {
    hr(); cat("WARNINGS:\n"); for (w in x$warnings) cat("  ! ", w, "\n", sep = "")
  }
  hr(); cat(sprintf("Done in %.1fs\n", x$elapsed_secs))
  invisible(x)
}
