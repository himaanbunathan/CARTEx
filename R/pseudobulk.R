# Pseudobulk aggregation, following the reference CARTEx implementation
# (vandontduong/CARTEx): cells are grouped (e.g. by cluster / sample), counts
# are summed within each group with Seurat::AggregateExpression-style behaviour,
# the pseudobulk profiles are LogNormalized, and then scored exactly like single
# cells (weighted dot product, Z-scored across the pseudobulk samples). An
# optional random partition into N pseudo-replicates per group mirrors the
# reference `PseudoBulkLabels()` helper.

#' Build pseudo-replicate bin labels (reference `PseudoBulkLabels`).
#' @keywords internal
#' @noRd
.pseudobulk_replicates <- function(n_cells, n_replicates, seed = 1L) {
  if (n_replicates <= 1) return(rep(1L, n_cells))
  withr_seed <- get0(".Random.seed", envir = globalenv())
  set.seed(seed)
  on.exit(if (!is.null(withr_seed)) assign(".Random.seed", withr_seed, envir = globalenv()))
  sample.int(n_replicates, n_cells, replace = TRUE)
}

#' Aggregate cells into pseudobulk profiles and score them.
#'
#' @param expr genes x cells matrix for the CARTEx genes (counts when
#'   `is_counts = TRUE`, otherwise log-normalized).
#' @param cell_total Per-cell library size over ALL genes (used to normalize
#'   summed counts); required when `is_counts = TRUE`.
#' @param group_keys Character vector (length = n cells) giving each cell's group.
#' @param weights Numeric CARTEx weight per gene (row-aligned to `expr`).
#' @param is_counts Whether `expr` holds raw counts (sum+normalize) or
#'   normalized values (average).
#' @param scale_factor LogNormalize scale factor.
#' @return list(scores = data.frame per pseudobulk sample, group = data.frame
#'   with n_cells, action = description).
#' @keywords internal
#' @noRd
.score_pseudobulk <- function(expr, cell_total, group_keys, weights,
                              is_counts = TRUE, scale_factor = 1e4) {
  group_keys <- as.character(group_keys)
  keep <- !is.na(group_keys)
  expr <- expr[, keep, drop = FALSE]
  group_keys <- group_keys[keep]
  if (!is.null(cell_total)) cell_total <- cell_total[keep]

  grp <- factor(group_keys)
  ind <- Matrix::fac2sparse(grp)                 # levels x cells indicator
  n_per <- as.numeric(Matrix::rowSums(ind))      # cells per group

  pb <- expr %*% Matrix::t(ind)                  # genes x groups (sum)
  pb <- as.matrix(pb)
  colnames(pb) <- levels(grp)

  if (is_counts) {
    if (is.null(cell_total)) stop("Pseudobulk of counts needs library sizes.")
    pb_lib <- as.numeric(ind %*% cell_total)     # summed library size per group
    pb_lib[pb_lib == 0] <- 1
    pb_norm <- log1p(sweep(pb, 2, pb_lib, "/") * scale_factor)
    action <- sprintf("summed counts per group -> LogNormalize(%g)", scale_factor)
  } else {
    pb_norm <- sweep(pb, 2, n_per, "/")          # mean normalized expression
    action <- "averaged log-normalized expression per group"
  }

  raw <- as.numeric(crossprod(pb_norm, as.numeric(weights)))
  names(raw) <- colnames(pb_norm)
  z <- .zscore(raw)
  scores <- data.frame(
    pseudobulk_group = colnames(pb_norm),
    cartex_score = z,
    cartex_score_raw = raw,
    cartex_score_int = .integerize(z),
    n_cells = n_per,
    row.names = colnames(pb_norm),
    stringsAsFactors = FALSE
  )
  group_df <- data.frame(pseudobulk_group = levels(grp), n_cells = n_per,
                         row.names = levels(grp), stringsAsFactors = FALSE)
  list(scores = scores, group = group_df, action = action,
       n_pseudobulk = nlevels(grp))
}

#' Integerize a score: round and cap at +/- 4 (reference `integerize`).
#' @keywords internal
#' @noRd
.integerize <- function(score) {
  s <- round(score)
  s[s <= -4] <- -4
  s[s >=  4] <-  4
  s
}
