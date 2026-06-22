# CD8+ T-cell identification.
#
# CARTEx is only biologically valid for CD8+ T cells, so before scoring we
# classify each cell. There is no single perfect signal, so we use both:
#   1. a curated cell-type annotation column when one exists (most reliable);
#   2. marker-gene gating as a fallback (CD8A/CD8B positive, CD3 positive,
#      CD4 negative) -- which must require T-cell identity, because NK cells
#      also express CD8A.
# Cells that cannot be resolved confidently (marker dropout, doublets, or a
# generic "T cell" label) are labelled "ambiguous_*" rather than forced into a
# call. Nothing is dropped here; the orchestrator decides whether to flag or
# restrict. Every decision is summarized for the report.

#' Marker genes used for CD8+ T-cell gating (symbol + Ensembl, for robustness).
#' @keywords internal
#' @noRd
.cd8_marker_map <- function() {
  data.frame(
    role   = c("CD8", "CD8", "CD4", "T", "T", "T"),
    symbol = c("CD8A", "CD8B", "CD4", "CD3D", "CD3E", "CD3G"),
    ensembl_id = c("ENSG00000153563", "ENSG00000172116", "ENSG00000010610",
                   "ENSG00000167286", "ENSG00000198851", "ENSG00000160654"),
    stringsAsFactors = FALSE
  )
}

#' Candidate metadata columns that may hold a cell-type annotation.
#' @keywords internal
#' @noRd
.celltype_column_candidates <- function() {
  c("cell_type", "celltype", "cell.type", "CellType", "cell_type_annotation",
    "annotation", "annotations", "seurat_annotations", "predicted.celltype",
    "predicted.id", "major_celltype", "cluster_annotation", "labels", "ident")
}

#' Classify cells as CD8+ T cells.
#'
#' @param meta Cell metadata data.frame.
#' @param marker_expr Optional markers x cells matrix (rows named by symbol).
#' @param method One of "auto", "metadata", "markers", "none".
#' @param cd8_column Optional explicit annotation column name.
#' @return list(status = factor per cell, method_used, column, counts, notes).
#' @keywords internal
#' @noRd
.classify_cd8 <- function(meta, marker_expr = NULL, method = "auto",
                          cd8_column = NULL) {
  n <- if (!is.null(marker_expr)) ncol(marker_expr) else nrow(meta)
  cell_names <- if (!is.null(marker_expr)) colnames(marker_expr) else rownames(meta)

  find_annotation_column <- function() {
    if (!is.null(cd8_column)) {
      if (!cd8_column %in% names(meta))
        stop("cd8_column '", cd8_column, "' not found in cell metadata.")
      return(cd8_column)
    }
    cand <- intersect(.celltype_column_candidates(), names(meta))
    # Prefer a column whose values actually mention CD8/CD4/T cells.
    for (col in cand) {
      vals <- as.character(meta[[col]])
      if (any(grepl("CD8|CD4|T cell|T-cell|T lymphocyte", vals, ignore.case = TRUE)))
        return(col)
    }
    if (length(cand)) cand[1] else NA_character_
  }

  classify_by_metadata <- function(col) {
    lab <- as.character(meta[[col]])
    status <- rep("unknown", length(lab))
    is_cd8 <- grepl("CD8", lab, ignore.case = TRUE)
    is_cd4 <- grepl("CD4", lab, ignore.case = TRUE)
    is_generic_T <- grepl("T cell|T-cell|T lymphocyte|\\bTcell\\b", lab, ignore.case = TRUE)
    is_nonT <- grepl(paste0("B cell|B-cell|plasma|NK|natural killer|innate lymphoid|",
                            "\\bILC\\b|myeloid|monocyte|macrophage|dendritic|\\bDC\\b|",
                            "neutrophil|granulocyte|mast|erythro|platelet|megakaryo|",
                            "epitheli|fibroblast|endotheli|stromal"),
                     lab, ignore.case = TRUE)
    status[is_nonT] <- "CD8_neg"
    status[is_cd4 & !is_cd8] <- "CD8_neg"
    status[is_generic_T & !is_cd8 & !is_cd4] <- "ambiguous_generic"
    status[is_cd8 & !is_cd4] <- "CD8_pos"
    status[is_cd8 & is_cd4] <- "ambiguous_doublet"
    status
  }

  classify_by_markers <- function() {
    get <- function(sym) {
      if (sym %in% rownames(marker_expr)) as.numeric(marker_expr[sym, ]) else rep(0, n)
    }
    cd8 <- (get("CD8A") > 0) | (get("CD8B") > 0)
    cd4 <- get("CD4") > 0
    tcell <- (get("CD3D") > 0) | (get("CD3E") > 0) | (get("CD3G") > 0)
    status <- rep("unknown", n)
    status[!tcell] <- "non_T"
    status[tcell & !cd8 & !cd4] <- "ambiguous_dropout"
    status[tcell & cd8 &  cd4]  <- "ambiguous_doublet"
    status[tcell & !cd8 &  cd4] <- "CD8_neg"
    status[tcell &  cd8 & !cd4] <- "CD8_pos"
    status
  }

  notes <- character(0)
  col <- NA_character_
  if (method == "none") {
    status <- rep("not_assessed", n)
    method_used <- "none"
  } else if (method == "metadata") {
    col <- find_annotation_column()
    if (is.na(col)) stop("No cell-type annotation column found; try cd8_method='markers'.")
    status <- classify_by_metadata(col); method_used <- "metadata"
  } else if (method == "markers") {
    if (is.null(marker_expr)) stop("Marker expression unavailable for marker-based CD8 gating.")
    status <- classify_by_markers(); method_used <- "markers"
  } else { # auto
    col <- find_annotation_column()
    if (!is.na(col)) {
      status <- classify_by_metadata(col); method_used <- "metadata"
      notes <- c(notes, sprintf("used annotation column '%s'", col))
    } else if (!is.null(marker_expr)) {
      status <- classify_by_markers(); method_used <- "markers"
      notes <- c(notes, "no annotation column; gated on marker genes")
    } else {
      status <- rep("unknown", n); method_used <- "none"
      notes <- c(notes, "no annotation column and no marker genes available")
    }
  }

  status <- factor(status, levels = c("CD8_pos", "CD8_neg", "ambiguous_generic",
                                      "ambiguous_dropout", "ambiguous_doublet",
                                      "non_T", "unknown", "not_assessed"))
  names(status) <- cell_names
  list(status = status, method_used = method_used, column = col,
       counts = table(status), notes = notes)
}
