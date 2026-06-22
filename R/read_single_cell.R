# Readers for single-cell objects in several on-disk formats.
#
# Goal: score huge files on a laptop. For the HDF5-backed formats
# (.h5ad / .loom / .h5Seurat) we extract ONLY the genes we need -- the CARTEx
# signature genes plus a handful of CD8 marker genes -- from the sparse matrix
# in a single streaming pass, never materializing the full genes x cells dense
# matrix or a whole Seurat object. Per-cell library sizes (needed to
# LogNormalize raw counts and to build pseudobulk profiles) are accumulated in
# the same pass, for free.
#
# HDF5 access uses rhdf5 (Bioconductor), which bundles its own HDF5 library and
# therefore installs cleanly on Mac/Windows/Linux without system dependencies.

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Marker genes shaped like the CARTEx gene map, so the same matcher works.
#' @keywords internal
#' @noRd
.marker_gene_map <- function() {
  m <- .cd8_marker_map()
  data.frame(cartex_symbol = m$symbol, clean_symbol = m$symbol,
             ensembl_id = m$ensembl_id, weight = NA_real_,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Low-level sparse streaming
# ---------------------------------------------------------------------------

#' Stream a cell-major sparse matrix, keeping only target gene rows.
#'
#' Works for AnnData CSR (`csr_matrix`, indptr per cell, indices = gene) and for
#' Seurat dgCMatrix data (CSC with columns = cells, indptr per cell, indices =
#' gene). Both are "cell-major with gene minor indices".
#'
#' @return list(expr = genes x cells dgCMatrix in `gene_index` order,
#'   cell_total = per-cell sum over ALL genes of this layer).
#' @keywords internal
#' @noRd
.sparse_extract_cellmajor <- function(path, data_path, indices_path, indptr_path,
                                      gene_index, n_genes, n_cells, chunk = 4000L) {
  indptr <- as.numeric(rhdf5::h5read(path, indptr_path))
  keep <- logical(n_genes); keep[gene_index] <- TRUE
  out_pos <- integer(n_genes); out_pos[gene_index] <- seq_along(gene_index)

  i_list <- list(); j_list <- list(); x_list <- list(); blk <- 0L
  cell_total <- numeric(n_cells)

  c0 <- 1L
  while (c0 <= n_cells) {
    c1 <- min(c0 + chunk - 1L, n_cells)
    a <- indptr[c0] + 1L; b <- indptr[c1 + 1L]
    if (b >= a) {
      d  <- as.numeric(rhdf5::h5read(path, data_path,    index = list(a:b)))
      gi <- as.integer(rhdf5::h5read(path, indices_path, index = list(a:b)))
      nper <- diff(indptr[c0:(c1 + 1L)])
      cell_ids <- rep.int(c0:c1, nper)
      cell_total[c0:c1] <- as.numeric(tapply(d, factor(cell_ids, levels = c0:c1),
                                             sum, default = 0))
      sel <- keep[gi + 1L]
      if (any(sel)) {
        blk <- blk + 1L
        i_list[[blk]] <- out_pos[gi[sel] + 1L]
        j_list[[blk]] <- cell_ids[sel]
        x_list[[blk]] <- d[sel]
      }
    }
    c0 <- c1 + 1L
  }
  cell_total[is.na(cell_total)] <- 0
  expr <- Matrix::sparseMatrix(
    i = if (blk) unlist(i_list) else integer(0),
    j = if (blk) unlist(j_list) else integer(0),
    x = if (blk) unlist(x_list) else numeric(0),
    dims = c(length(gene_index), n_cells))
  list(expr = expr, cell_total = cell_total)
}

#' Extract specific gene rows from an AnnData CSC matrix (gene-major).
#' @keywords internal
#' @noRd
.sparse_extract_genemajor <- function(path, data_path, indices_path, indptr_path,
                                      gene_index, n_cells, need_total) {
  indptr <- as.numeric(rhdf5::h5read(path, indptr_path))
  i_list <- list(); j_list <- list(); x_list <- list()
  for (r in seq_along(gene_index)) {
    g <- gene_index[r]
    a <- indptr[g] + 1L; b <- indptr[g + 1L]
    if (b >= a) {
      d  <- as.numeric(rhdf5::h5read(path, data_path,    index = list(a:b)))
      ci <- as.integer(rhdf5::h5read(path, indices_path, index = list(a:b))) + 1L
      i_list[[r]] <- rep.int(r, length(d)); j_list[[r]] <- ci; x_list[[r]] <- d
    }
  }
  expr <- Matrix::sparseMatrix(i = unlist(i_list), j = unlist(j_list), x = unlist(x_list),
                               dims = c(length(gene_index), n_cells))
  cell_total <- NULL
  if (need_total) {
    cell_total <- numeric(n_cells)
    n_genes <- length(indptr) - 1L
    for (g in seq_len(n_genes)) {
      a <- indptr[g] + 1L; b <- indptr[g + 1L]
      if (b >= a) {
        d  <- as.numeric(rhdf5::h5read(path, data_path,    index = list(a:b)))
        ci <- as.integer(rhdf5::h5read(path, indices_path, index = list(a:b))) + 1L
        cell_total[ci] <- cell_total[ci] + d
      }
    }
  }
  list(expr = expr, cell_total = cell_total)
}

# ---------------------------------------------------------------------------
# AnnData (.h5ad)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.h5ad_features <- function(path) {
  ids <- as.character(rhdf5::h5read(path, "var/_index"))
  fn <- tryCatch(rhdf5::h5read(path, "var/feature_name"), error = function(e) NULL)
  if (is.list(fn) && all(c("categories", "codes") %in% names(fn))) {
    symbols <- as.character(fn$categories)[fn$codes + 1L]
  } else if (!is.null(fn)) {
    symbols <- as.character(fn)
  } else {
    symbols <- ids
  }
  list(symbols = symbols, ids = ids)
}

#' @keywords internal
#' @noRd
.h5ad_read_dataframe <- function(path, group) {
  raw <- rhdf5::h5read(path, group)
  idx <- raw[["_index"]]
  out <- list()
  for (nm in setdiff(names(raw), c("_index", "__categories"))) {
    v <- raw[[nm]]
    if (is.list(v) && all(c("categories", "codes") %in% names(v))) {
      out[[nm]] <- as.character(v$categories)[v$codes + 1L]
    } else if (is.list(v)) {
      next
    } else if (length(v) == length(idx)) {
      out[[nm]] <- as.vector(v)
    }
  }
  if (!length(out)) return(data.frame(row.names = as.character(idx)))
  data.frame(out, row.names = as.character(idx),
             stringsAsFactors = FALSE, check.names = FALSE)
}

#' Read an existing UMAP embedding (obsm/X_umap) from an h5ad, if present.
#' @keywords internal
#' @noRd
.h5ad_umap <- function(path, cells) {
  obsm <- tryCatch(rhdf5::h5ls(path, recursive = FALSE), error = function(e) NULL)
  has_obsm <- !is.null(obsm) && "obsm" %in% obsm$name
  if (!has_obsm) return(NULL)
  keys <- tryCatch(rhdf5::h5ls(path)$name[rhdf5::h5ls(path)$group == "/obsm"],
                   error = function(e) character(0))
  key <- intersect(c("X_umap", "X_UMAP", "umap", "UMAP"), keys)
  if (!length(key)) return(NULL)
  m <- tryCatch(rhdf5::h5read(path, paste0("obsm/", key[1])), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  m <- as.matrix(m)
  if (nrow(m) != length(cells) && ncol(m) == length(cells)) m <- t(m)
  if (nrow(m) != length(cells)) return(NULL)
  df <- data.frame(UMAP_1 = m[, 1], UMAP_2 = m[, 2], row.names = cells)
  df
}

#' Read specific gene rows from the best-matching h5ad layer.
#' @param want "normalized" (prefer log-normalized) or "counts" (prefer raw).
#' @keywords internal
#' @noRd
.h5ad_extract <- function(path, gene_index, want = "normalized", layer = NULL,
                          sample_n = 2e5L) {
  top <- rhdf5::h5ls(path, recursive = FALSE)$name
  shape <- rhdf5::h5readAttributes(path, "X")[["shape"]]
  if (is.null(shape)) {
    n_var <- length(rhdf5::h5read(path, "var/_index"))
    n_cells <- length(rhdf5::h5read(path, "X/indptr")) - 1L
    n_genes <- n_var
  } else {
    n_cells <- as.integer(shape[1]); n_genes <- as.integer(shape[2])
  }
  declared <- tryCatch(as.character(rhdf5::h5read(path, "uns/X_name")),
                       error = function(e) NULL)

  candidates <- if (!is.null(layer)) layer else c("X", if ("raw" %in% top) "raw/X")
  inspected <- lapply(candidates, function(lp) {
    s <- tryCatch(as.numeric(rhdf5::h5read(path, paste0(lp, "/data"),
                                           index = list(1:min(sample_n, 5e5L)))),
                  error = function(e) NULL)
    if (is.null(s)) return(NULL)
    info <- .inspect_expression(s, declared = if (lp == "X") declared else lp)
    c(list(layer = lp), info)
  })
  inspected <- Filter(Negate(is.null), inspected)
  if (!length(inspected)) stop("Could not read any expression layer from ", path)

  pref <- if (want == "counts") c("counts", "log_normalized", "linear_normalized")
          else c("log_normalized", "linear_normalized", "counts")
  pick <- NULL
  for (st in pref) {
    hit <- Filter(function(x) x$state == st, inspected)
    if (length(hit)) { pick <- hit[[1]]; break }
  }
  if (is.null(pick)) pick <- inspected[[1]]

  lp <- pick$layer
  enc <- rhdf5::h5readAttributes(path, lp)[["encoding-type"]] %||% "csr_matrix"
  need_total <- TRUE   # always compute library sizes (cheap; needed downstream)
  if (enc == "csc_matrix") {
    ext <- .sparse_extract_genemajor(path, paste0(lp, "/data"), paste0(lp, "/indices"),
                                     paste0(lp, "/indptr"), gene_index, n_cells, need_total)
  } else {
    ext <- .sparse_extract_cellmajor(path, paste0(lp, "/data"), paste0(lp, "/indices"),
                                     paste0(lp, "/indptr"), gene_index, n_genes, n_cells)
  }
  list(expr = ext$expr, cell_total = ext$cell_total, state = pick$state,
       declared = declared, layer_used = lp, inspect_reason = pick$reason,
       n_cells = n_cells, n_genes = n_genes)
}

# ---------------------------------------------------------------------------
# Helpers to split a unioned extraction back into CARTEx vs marker matrices
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.split_targets <- function(expr_union, union_index, cartex_match, marker_match) {
  cm <- cartex_match$matched
  cartex_expr <- expr_union[match(cm$feature_index, union_index), , drop = FALSE]
  rownames(cartex_expr) <- cm$cartex_symbol
  marker_expr <- NULL
  if (!is.null(marker_match) && marker_match$n_matched > 0) {
    mk <- marker_match$matched
    marker_expr <- expr_union[match(mk$feature_index, union_index), , drop = FALSE]
    rownames(marker_expr) <- mk$cartex_symbol
  }
  list(cartex_expr = cartex_expr, marker_expr = marker_expr)
}

# ---------------------------------------------------------------------------
# Loom (.loom)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.loom_load <- function(path, gene_map, marker_map, want, layer = NULL) {
  ra <- rhdf5::h5read(path, "row_attrs")
  ca <- rhdf5::h5read(path, "col_attrs")
  symbols <- as.character(ra[["Gene"]] %||% ra[["gene_names"]] %||% ra[[1]])
  ids <- if (!is.null(ra[["Accession"]])) as.character(ra[["Accession"]]) else NULL
  cm <- .match_cartex_genes(symbols, ids, gene_map)
  mk <- .match_cartex_genes(symbols, ids, marker_map)
  union_index <- sort(unique(c(cm$matched$feature_index, mk$matched$feature_index)))

  mat_path <- if (!is.null(layer)) paste0("layers/", layer) else "matrix"
  sub <- rhdf5::h5read(path, mat_path, index = list(union_index, NULL))
  expr_union <- Matrix::Matrix(sub, sparse = TRUE)
  cells <- as.character(ca[["CellID"]] %||% ca[["obs_names"]] %||% seq_len(ncol(expr_union)))
  colnames(expr_union) <- cells
  sp <- .split_targets(expr_union, union_index, cm, mk)

  meta_cols <- Filter(function(v) is.atomic(v) && length(v) == ncol(expr_union),
                      ca[setdiff(names(ca), "CellID")])
  meta <- if (length(meta_cols))
    data.frame(meta_cols, row.names = cells, stringsAsFactors = FALSE, check.names = FALSE)
  else data.frame(row.names = cells)

  st <- .inspect_expression(as.numeric(sp$cartex_expr@x))$state
  list(cartex_expr = sp$cartex_expr, marker_expr = sp$marker_expr, matched = cm,
       cell_total = Matrix::colSums(expr_union), state = st, declared = NULL,
       layer_used = mat_path, meta = meta, n_cells = ncol(expr_union),
       n_genes = length(symbols))
}

# ---------------------------------------------------------------------------
# h5Seurat (.h5Seurat)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.h5seurat_load <- function(path, gene_map, marker_map, want, assay = NULL, layer = NULL) {
  top <- rhdf5::h5ls(path, recursive = FALSE)
  if (!"assays" %in% top$name)
    stop(".h5Seurat file has no 'assays' group; convert with SeuratDisk first.")
  assay_names <- rhdf5::h5read(path, "assays")
  if (is.null(assay)) assay <- names(assay_names)[1]
  base <- paste0("assays/", assay)
  feats <- as.character(rhdf5::h5read(path, paste0(base, "/features")))
  slot <- layer %||% if (want == "counts") "counts" else "data"
  cm <- .match_cartex_genes(feats, NULL, gene_map)
  mk <- .match_cartex_genes(feats, NULL, marker_map)
  union_index <- sort(unique(c(cm$matched$feature_index, mk$matched$feature_index)))
  dpath <- paste0(base, "/", slot)
  n_genes <- length(feats)
  n_cells <- length(as.numeric(rhdf5::h5read(path, paste0(dpath, "/indptr")))) - 1L
  ext <- .sparse_extract_cellmajor(path, paste0(dpath, "/data"),
                                   paste0(dpath, "/indices"), paste0(dpath, "/indptr"),
                                   union_index, n_genes, n_cells)
  cells <- tryCatch(as.character(rhdf5::h5read(path, "cell.names")),
                    error = function(e) as.character(seq_len(n_cells)))
  colnames(ext$expr) <- cells
  sp <- .split_targets(ext$expr, union_index, cm, mk)
  meta <- tryCatch(.h5ad_read_dataframe(path, "meta.data"),
                   error = function(e) data.frame(row.names = cells))
  list(cartex_expr = sp$cartex_expr, marker_expr = sp$marker_expr, matched = cm,
       cell_total = ext$cell_total,
       state = if (slot == "counts") "counts" else "log_normalized",
       declared = slot, layer_used = paste0(assay, "/", slot), meta = meta,
       n_cells = n_cells, n_genes = n_genes)
}

# ---------------------------------------------------------------------------
# In-memory Seurat object / plain matrix
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.subset_rows <- function(mat, match_res) {
  m <- match_res$matched
  if (!nrow(m)) return(NULL)
  out <- mat[m$feature_index, , drop = FALSE]
  rownames(out) <- m$cartex_symbol
  out
}

#' Pull an existing UMAP embedding from a Seurat object, if present.
#' @keywords internal
#' @noRd
.seurat_umap <- function(obj) {
  reds <- tryCatch(names(obj@reductions), error = function(e) NULL)
  key <- intersect(c("umap", "UMAP", "X_umap", "umap.cca", "ref.umap"), reds)
  if (length(key)) {
    emb <- SeuratObject::Embeddings(obj, reduction = key[1])
    if (ncol(emb) >= 2)
      return(data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                        row.names = rownames(emb)))
  }
  # Fall back to UMAP_1/UMAP_2-style metadata columns.
  md <- obj@meta.data
  cols <- grep("^(umap|UMAP)[_.]?1$", names(md), value = TRUE)
  if (length(cols)) {
    c2 <- sub("1$", "2", cols[1])
    if (c2 %in% names(md))
      return(data.frame(UMAP_1 = md[[cols[1]]], UMAP_2 = md[[c2]],
                        row.names = rownames(md)))
  }
  NULL
}

#' @keywords internal
#' @noRd
.from_seurat <- function(obj, gene_map, marker_map, want, assay = NULL, layer = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("Reading Seurat objects requires the 'SeuratObject' package.")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(obj)
  layers_avail <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  slot <- layer %||% if (want == "counts" && (is.null(layers_avail) || "counts" %in% layers_avail))
    "counts" else if (!is.null(layers_avail) && "data" %in% layers_avail) "data" else "counts"
  mat <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = slot),
                  error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = slot))
  feats <- rownames(mat)
  cm <- .match_cartex_genes(feats, NULL, gene_map)
  mk <- .match_cartex_genes(feats, NULL, marker_map)
  state <- if (slot == "counts") "counts" else "log_normalized"
  list(cartex_expr = .subset_rows(mat, cm), marker_expr = .subset_rows(mat, mk),
       matched = cm, cell_total = Matrix::colSums(mat), state = state,
       declared = slot, layer_used = paste0(assay, "/", slot),
       meta = obj@meta.data, umap = .seurat_umap(obj),
       n_cells = ncol(mat), n_genes = nrow(mat), seurat = obj)
}

#' @keywords internal
#' @noRd
.from_matrix <- function(mat, gene_map, marker_map) {
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  feats <- rownames(mat)
  if (is.null(feats)) stop("Matrix input must have gene symbols as row names.")
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  cm <- .match_cartex_genes(feats, NULL, gene_map)
  mk <- .match_cartex_genes(feats, NULL, marker_map)
  cartex_expr <- .subset_rows(mat, cm)
  st <- if (is.null(cartex_expr)) "unknown"
        else .inspect_expression(as.numeric(cartex_expr@x))$state
  list(cartex_expr = cartex_expr, marker_expr = .subset_rows(mat, mk),
       matched = cm, cell_total = Matrix::colSums(mat), state = st,
       declared = NULL, layer_used = "matrix",
       meta = data.frame(row.names = colnames(mat)),
       n_cells = ncol(mat), n_genes = nrow(mat))
}

# ---------------------------------------------------------------------------
# Dispatcher used by cartex_sc()
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.load_expression <- function(input, gene_map, want = "normalized",
                             assay = NULL, layer = NULL) {
  marker_map <- .marker_gene_map()

  finish <- function(res, fmt, src) {
    res$provenance <- list(format = fmt, source = src,
                           assay = assay %||% NA, layer = res$layer_used)
    res
  }

  if (inherits(input, "Seurat"))
    return(finish(.from_seurat(input, gene_map, marker_map, want, assay, layer),
                  "Seurat (in-memory)", "object"))
  if (is.matrix(input) || is.data.frame(input) || inherits(input, "Matrix"))
    return(finish(.from_matrix(input, gene_map, marker_map), "matrix (in-memory)", "object"))
  if (!is.character(input) || length(input) != 1)
    stop("`input` must be a file path, a Seurat object, or a gene x cell matrix.")
  if (!file.exists(input)) stop("File not found: ", input)

  ext <- tolower(tools::file_ext(input))
  if (ext == "rds") {
    obj <- readRDS(input)
    res <- if (inherits(obj, "Seurat")) .from_seurat(obj, gene_map, marker_map, want, assay, layer)
           else .from_matrix(obj, gene_map, marker_map)
  } else if (ext %in% c("rdata", "rda")) {
    e <- new.env(); load(input, envir = e); obj <- get(ls(e)[1], envir = e)
    res <- if (inherits(obj, "Seurat")) .from_seurat(obj, gene_map, marker_map, want, assay, layer)
           else .from_matrix(obj, gene_map, marker_map)
  } else if (ext == "h5ad") {
    .require_rhdf5()
    feats <- .h5ad_features(input)
    cm <- .match_cartex_genes(feats$symbols, feats$ids, gene_map)
    mk <- .match_cartex_genes(feats$symbols, feats$ids, marker_map)
    union_index <- sort(unique(c(cm$matched$feature_index, mk$matched$feature_index)))
    ex <- .h5ad_extract(input, union_index, want = want, layer = layer)
    colnames(ex$expr) <- as.character(rhdf5::h5read(input, "obs/_index"))
    sp <- .split_targets(ex$expr, union_index, cm, mk)
    meta <- .h5ad_read_dataframe(input, "obs")
    umap <- .h5ad_umap(input, colnames(sp$cartex_expr))
    rhdf5::h5closeAll()
    res <- list(cartex_expr = sp$cartex_expr, marker_expr = sp$marker_expr, matched = cm,
                cell_total = ex$cell_total, state = ex$state, declared = ex$declared,
                layer_used = ex$layer_used, inspect_reason = ex$inspect_reason,
                meta = meta, umap = umap, n_cells = ex$n_cells, n_genes = ex$n_genes)
  } else if (ext == "loom") {
    .require_rhdf5(); res <- .loom_load(input, gene_map, marker_map, want, layer); rhdf5::h5closeAll()
  } else if (ext == "h5seurat") {
    .require_rhdf5()
    res <- tryCatch(.h5seurat_load(input, gene_map, marker_map, want, assay, layer),
                    error = function(e) { rhdf5::h5closeAll(); stop(conditionMessage(e)) })
    rhdf5::h5closeAll()
  } else {
    stop("Unsupported file extension: '", ext,
         "'. Supported: .rds, .RData, .h5ad, .h5Seurat, .loom.")
  }
  finish(res, ext, input)
}

#' @keywords internal
#' @noRd
.require_rhdf5 <- function() {
  if (!requireNamespace("rhdf5", quietly = TRUE))
    stop("Reading HDF5-backed formats (.h5ad/.loom/.h5Seurat) requires the ",
         "Bioconductor package 'rhdf5'. Install with:\n",
         "  install.packages('BiocManager'); BiocManager::install('rhdf5')",
         call. = FALSE)
}

#' Read single-cell expression for the CARTEx genes from any supported format
#'
#' A thin, user-facing wrapper around the internal format readers. It returns a
#' lightweight list (matched CARTEx genes x cells sparse matrix, CD8 marker
#' genes x cells, cell metadata, and provenance) without loading the entire
#' object into memory. Most users should call [cartex_sc()] instead, which also
#' scores the cells.
#'
#' @param path Path to a `.rds`, `.RData`, `.h5ad`, `.h5Seurat`, or `.loom`
#'   file (or an in-memory Seurat object / gene x cell matrix).
#' @param want Either `"normalized"` (default) or `"counts"`; which layer to
#'   prefer when both are available.
#' @param assay Optional assay name (Seurat / h5Seurat).
#' @param layer Optional layer/slot name to read instead of the auto-selected one.
#' @return A list with `cartex_expr`, `marker_expr`, `meta`, `matched`, `state`,
#'   and `provenance`.
#' @export
read_single_cell <- function(path, want = "normalized", assay = NULL, layer = NULL) {
  .load_expression(path, .cartex_gene_map(), want = want, assay = assay, layer = layer)
}
