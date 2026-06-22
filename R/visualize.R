# Visualization and HTML reporting for CARTEx single-cell results.
#
# cartex_umap()        -- side-by-side UMAP: cells coloured by annotation and by
#                         CARTEx score. Uses an existing UMAP embedding if the
#                         object already has one; otherwise (for a Seurat object)
#                         runs the standard Seurat pipeline with default
#                         parameters to compute one.
# cartex_html_report() -- a self-contained HTML summary of a cartex_sc() run,
#                         optionally embedding the UMAP figure.

#' @keywords internal
#' @noRd
.need_pkg <- function(pkg, what) {
  miss <- pkg[!vapply(pkg, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop(what, " requires package(s): ", paste(miss, collapse = ", "),
         ". Install them and retry.", call. = FALSE)
}

#' Compute a default UMAP for a Seurat object (standard pipeline).
#' @keywords internal
#' @noRd
.compute_umap_seurat <- function(obj, assay = NULL) {
  .need_pkg("Seurat", "Computing a UMAP")
  message("No UMAP embedding found; computing one with default Seurat parameters",
          " (NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA -> RunUMAP)...")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(obj)
  ld <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  if (is.null(ld) || !any(grepl("^data", ld)))
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = 1:30, verbose = FALSE)
  emb <- SeuratObject::Embeddings(obj, "umap")
  data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], row.names = rownames(emb))
}

#' Assemble a plotting data.frame (UMAP coords + score + annotation) from a
#' cartex_sc() result.
#' @keywords internal
#' @noRd
.build_plot_df <- function(x, input = NULL, annotation = NULL,
                           score = "cartex_score", compute = TRUE) {
  pick_annotation <- function(nms, md) {
    if (!is.null(annotation)) {
      if (!annotation %in% nms) stop("annotation column '", annotation, "' not found.")
      return(annotation)
    }
    rep <- attr(x, "cartex_report")
    if (!is.null(rep) && !is.null(rep$cd8$column) && !is.na(rep$cd8$column) &&
        rep$cd8$column %in% nms) return(rep$cd8$column)
    cand <- intersect(.celltype_column_candidates(), nms)
    if (length(cand)) cand[1] else "cd8_status"
  }

  grab <- function(src, col) if (col %in% names(src)) src[[col]] else NA

  if (inherits(x, "Seurat")) {
    md <- x@meta.data
    umap <- .seurat_umap(x)
    if (is.null(umap)) {
      if (!compute) stop("No UMAP in object and compute = FALSE.")
      umap <- .compute_umap_seurat(x)
    }
    ann_col <- pick_annotation(names(md), md)
    df <- data.frame(
      UMAP_1 = umap[rownames(md), "UMAP_1"], UMAP_2 = umap[rownames(md), "UMAP_2"],
      score = md[[score]], annotation = as.character(md[[ann_col]]),
      cd8_status = as.character(grab(md, "cd8_status")),
      raw = grab(md, "cartex_score_raw"),
      row.names = rownames(md), stringsAsFactors = FALSE)
  } else if (is.data.frame(x)) {
    nms <- names(x)
    if (all(c("UMAP_1", "UMAP_2") %in% nms)) {
      umap <- x[, c("UMAP_1", "UMAP_2")]
    } else if (!is.null(input) && is.character(input) &&
               tolower(tools::file_ext(input)) == "h5ad") {
      .require_rhdf5(); umap <- .h5ad_umap(input, rownames(x)); rhdf5::h5closeAll()
      if (is.null(umap)) stop("No UMAP_1/UMAP_2 columns and none found in '", input, "'.")
    } else {
      stop("No UMAP coordinates available. Provide a scored Seurat object, a ",
           "data.frame with UMAP_1/UMAP_2 columns, or `input` pointing to an ",
           ".h5ad with obsm/X_umap.")
    }
    if (!score %in% nms) stop("score column '", score, "' not found in the result.")
    ann_col <- pick_annotation(nms, x)
    ann <- if (ann_col %in% nms) as.character(x[[ann_col]]) else NA_character_
    df <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], score = x[[score]],
                     annotation = ann,
                     cd8_status = as.character(grab(x, "cd8_status")),
                     raw = grab(x, "cartex_score_raw"),
                     row.names = rownames(x), stringsAsFactors = FALSE)
  } else {
    stop("`x` must be a cartex_sc() result (Seurat object or data.frame).")
  }
  attr(df, "annotation_col") <- ann_col
  attr(df, "score_col") <- score
  df
}

#' Side-by-side UMAP: cell annotation and CARTEx score
#'
#' Produces a two-panel UMAP from a [cartex_sc()] result -- the left panel
#' coloured by cell-type annotation, the right by CARTEx score -- on identical
#' coordinates. If the object already contains a UMAP embedding (Seurat
#' reduction, AnnData `obsm/X_umap`, or `UMAP_1`/`UMAP_2` columns) it is reused;
#' otherwise, for a Seurat object, a UMAP is computed with default Seurat
#' parameters.
#'
#' @param x A result from [cartex_sc()]: a scored Seurat object, or the scored
#'   data.frame (single-cell level).
#' @param input Optional original input path (e.g. the `.h5ad`); used to fetch an
#'   existing UMAP when `x` is a data.frame without UMAP columns.
#' @param annotation Metadata column to colour the left panel by. Defaults to the
#'   CD8 annotation column from the report, else a detected cell-type column.
#' @param score Score column for the right panel (default `"cartex_score"`).
#' @param compute If `TRUE` (default), compute a UMAP for a Seurat object that
#'   lacks one. Has no effect for data.frame input (the full matrix is not loaded).
#' @param cd8_only If `TRUE` (default), the CARTEx score is shown on the right
#'   panel only for CD8+ T cells; every other cell type (CD8-negative, ambiguous,
#'   non-T, unknown) is greyed, because the CARTEx score is only biologically
#'   valid for CD8+ T cells and colouring other cells by it is misleading. All
#'   cells are still drawn on both panels, so the population context is preserved.
#'   Set `cd8_only = FALSE` to colour every cell by its score. Requires a
#'   `cd8_status` column (present in any `cartex_sc()` result).
#' @param renormalize If `TRUE` and `cd8_only = TRUE`, re-Z-score the raw
#'   per-cell score *within the CD8+ population* before plotting, so the colour
#'   scale reflects variation among the cells CARTEx is valid for rather than the
#'   whole-object distribution (which is diluted by unscored non-CD8 cells).
#'   Default `FALSE` (show the score column as computed, just masked).
#' @param na_colour Colour for cells that are not scored on the right panel
#'   (default "grey85").
#' @param point_size Point size (default 0.4).
#' @param file Optional path to save the figure (png/pdf via `ggplot2::ggsave`).
#' @return A patchwork/ggplot object (invisibly if `file` is given).
#' @export
cartex_umap <- function(x, input = NULL, annotation = NULL,
                        score = "cartex_score", compute = TRUE,
                        cd8_only = TRUE, renormalize = FALSE,
                        na_colour = "grey85",
                        point_size = 0.4, file = NULL) {
  .need_pkg(c("ggplot2", "patchwork"), "cartex_umap")
  df <- .build_plot_df(x, input, annotation, score, compute)
  ann_col <- attr(df, "annotation_col")

  score_title <- "CARTEx score"
  n_masked <- 0L
  if (cd8_only) {
    if (all(is.na(df$cd8_status)))
      stop("cd8_only = TRUE needs a 'cd8_status' column; none found in the ",
           "result. Pass cd8_only = FALSE to colour every cell by its score.")
    is_cd8 <- !is.na(df$cd8_status) & df$cd8_status == "CD8_pos"
    n_masked <- sum(!is_cd8)
    if (renormalize && !all(is.na(df$raw))) {
      # Re-Z-score the raw per-cell score within CD8+ cells only, so the colour
      # scale reflects variation among the cells CARTEx is valid for.
      z <- rep(NA_real_, nrow(df))
      z[is_cd8] <- .zscore(df$raw[is_cd8])
      df$score <- z
      score_title <- "CARTEx score\n(CD8+, re-normalized)"
    } else {
      # Mask: keep the computed score, but show it for CD8+ T cells only.
      df$score[!is_cd8] <- NA_real_
      score_title <- "CARTEx score\n(CD8+ T cells)"
    }
    if (n_masked)
      message(sprintf(
        "cartex_umap: %d of %d cells are not CD8+ T cells; greying them on the ",
        n_masked, nrow(df)),
        "score panel (CARTEx is only valid for CD8+ T cells). ",
        "Use cd8_only = FALSE to colour all cells.")
    # Draw greyed (non-CD8) cells first so scored CD8 cells sit on top.
    df <- df[order(!is.na(df$score)), ]
  }

  base <- ggplot2::ggplot(df, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  p_ann <- base +
    ggplot2::geom_point(ggplot2::aes(colour = .data$annotation), size = point_size) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 2.5), title = ann_col)) +
    ggplot2::labs(title = "Cell annotation")

  p_score <- base +
    ggplot2::geom_point(ggplot2::aes(colour = .data$score), size = point_size) +
    ggplot2::scale_colour_viridis_c(option = "C", name = score_title,
                                    na.value = na_colour) +
    ggplot2::labs(title = if (cd8_only) "CARTEx score (CD8+ only)" else "CARTEx score")

  combined <- patchwork::wrap_plots(p_ann, p_score, nrow = 1)
  if (!is.null(file)) {
    ggplot2::ggsave(file, combined, width = 12, height = 5, dpi = 150)
    message("Saved UMAP figure to ", file)
    return(invisible(combined))
  }
  combined
}

#' @keywords internal
#' @noRd
.plot_data_uri <- function(plot, width = 1100, height = 460) {
  .need_pkg(c("ggplot2", "base64enc"), "embedding a figure")
  tmp <- tempfile(fileext = ".png")
  ggplot2::ggsave(tmp, plot, width = width / 100, height = height / 100, dpi = 100)
  on.exit(unlink(tmp))
  base64enc::dataURI(file = tmp, mime = "image/png")
}

#' Write a self-contained HTML report for a CARTEx run
#'
#' Renders the [cartex_sc()] transparency report -- input provenance, gene
#' matching, CD8 composition, normalization decision, QC fields, score summary,
#' and any warnings -- to a standalone HTML file, optionally embedding the
#' side-by-side UMAP figure. No pandoc/rmarkdown required.
#'
#' @param x A result from [cartex_sc()] (with a `cartex_report` attribute).
#' @param file Output HTML path (default `"cartex_report.html"`).
#' @param umap Whether to embed the [cartex_umap()] figure (default TRUE).
#' @param input Optional original input path (for fetching an existing UMAP).
#' @param annotation Optional annotation column for the UMAP.
#' @param title Report title.
#' @param ... Further arguments forwarded to [cartex_umap()] (e.g.
#'   `cd8_only = TRUE`).
#' @return The output path, invisibly.
#' @export
cartex_html_report <- function(x, file = "cartex_report.html", umap = TRUE,
                               input = NULL, annotation = NULL,
                               title = "CARTEx single-cell report", ...) {
  .need_pkg("htmltools", "cartex_html_report")
  rep <- attr(x, "cartex_report")
  if (is.null(rep)) stop("`x` has no 'cartex_report' attribute; pass a cartex_sc() result.")
  tg <- htmltools::tags

  kv <- function(...) {
    rows <- list(...)
    tg$table(class = "kv",
      lapply(rows, function(r) tg$tr(tg$td(tg$b(r[[1]])), tg$td(r[[2]]))))
  }
  gm <- rep$gene_matching
  by_rule <- paste(sprintf("%s=%d", names(gm$by_rule), gm$by_rule), collapse = ", ")
  cd <- rep$cd8
  comp <- cd$counts[cd$counts > 0]
  comp_str <- paste(sprintf("%s=%d", names(comp), comp), collapse = ", ")

  sections <- list(
    tg$h1(title),
    tg$p(class = "muted", format(Sys.time())),
    tg$h2("Run summary"),
    kv(list("Input", rep$provenance$format),
       list("Source", as.character(rep$provenance$source %||% "")),
       list("Analysis level", rep$level),
       list("Cells in object", format(rep$dims$n_cells, big.mark = ",")),
       list("Units scored", sprintf("%s %s", format(rep$dims$n_scored, big.mark = ","),
                                    if (rep$level == "pseudobulk") "pseudobulk samples" else "cells")),
       list("Layer used", rep$normalization$layer)),
    tg$h2("Gene matching"),
    kv(list("Matched", sprintf("%d / %d (%.0f%%)", gm$n_matched, gm$n_requested,
                               100 * gm$match_fraction)),
       list("By rule", by_rule),
       list("Unmatched", paste(gm$unmatched, collapse = ", "))),
    tg$h2("CD8+ T-cell classification"),
    kv(list("Method", cd$method),
       list("Annotation column", as.character(cd$column %||% "n/a")),
       list("Mode", cd$mode),
       list("Composition", comp_str)),
    tg$h2("Normalization"),
    kv(list("Detected", rep$normalization$detected),
       list("Action", rep$normalization$action),
       list("Basis", as.character(rep$normalization$reason %||% ""))),
    tg$h2("QC metadata"),
    if (length(rep$qc$detected))
      tg$table(class = "kv",
        lapply(names(rep$qc$detected), function(k) {
          d <- rep$qc$detected[[k]]
          tg$tr(tg$td(tg$b(k)), tg$td(sprintf("[%s] min=%.3g median=%.3g max=%.3g",
                                              d$column, d$min, d$median, d$max)))
        }))
    else tg$p(class = "muted",
              "No nFeature/nCount/percent.mt columns found -- confirm QC upstream."),
    tg$h2("Score summary (raw weighted score)"),
    tg$pre(paste(utils::capture.output(print(rep$score_summary)), collapse = "\n"))
  )

  if (!is.null(rep$pseudobulk)) {
    pb <- rep$pseudobulk
    sections <- c(sections, list(
      tg$h2("Pseudobulk"),
      kv(list("Grouped by", paste(pb$columns, collapse = " + ")),
         list("Source", pb$source),
         list("Samples", pb$n_pseudobulk),
         list("Replicates", pb$replicates))))
  }

  if (isTRUE(umap)) {
    fig <- tryCatch({
      pl <- cartex_umap(x, input = input, annotation = annotation, ...)
      tg$div(tg$h2("UMAP"),
             tg$img(src = .plot_data_uri(pl), style = "max-width:100%;"))
    }, error = function(e) tg$p(class = "warn",
        paste("UMAP not rendered:", conditionMessage(e))))
    sections <- c(sections, list(fig))
  }

  sections <- c(sections, list(
    tg$h2("Caveat"), tg$p(rep$caveat),
    if (length(rep$warnings))
      tg$div(tg$h2("Warnings"),
             tg$ul(lapply(rep$warnings, function(w) tg$li(w))))
  ))

  css <- "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;
    margin:2rem auto;max-width:900px;color:#222;line-height:1.45}
    h1{margin-bottom:0}h2{border-bottom:1px solid #ddd;padding-bottom:4px;margin-top:1.6rem}
    table.kv{border-collapse:collapse}table.kv td{padding:2px 14px 2px 0;vertical-align:top}
    .muted{color:#888}.warn{color:#b00}pre{background:#f6f8fa;padding:10px;border-radius:6px}"

  page <- tg$html(tg$head(tg$meta(charset = "utf-8"), tg$title(title),
                          tg$style(htmltools::HTML(css))),
                  tg$body(sections))
  htmltools::save_html(page, file = file)
  message("Wrote HTML report to ", file)
  invisible(file)
}
