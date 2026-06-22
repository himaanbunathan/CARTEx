#' Calculate cartex score described in the paper
#'
#'
#' @param gene_counts A normalized gene counts data frame with genes as rows and samples as columns.
#' @param metadata An optional metadata dataframe corresponding to samples in the gene count dataframe
#' @return A dataframe of cartex scores.
#' @export
#'
cartex_score <- function(gene_counts, metadata=NULL) {

  if (!is.data.frame(gene_counts)) {
    stop("Input gene_counts must be a dataframe.")
  }

  if (is.null(rownames(gene_counts))) {
    stop("The rownames of Input dataframe should be genes.")
  }

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      stop("Input metadata must be a dataframe.")
    } else {
      metadata=metadata[order(rownames(metadata)),]
      gene_counts=gene_counts[,order(colnames(gene_counts))]
    }
  }

  precomputed_weights <- system.file("extdata", "weights.rds", package = "CARTEx")
  if (precomputed_weights == "") {
    stop("precomputed weights file not found.")
  }

  cartex_200 <- readRDS(precomputed_weights)

  common <- intersect(names(cartex_200), rownames(gene_counts))
  print(paste0("Num of genes used for cartex score = ", length(common)))

  Zscore <- function(s) {
    s <- as.numeric(s)
    z <- (s - mean(s)) / sd(s)
    return(z)
  }

  if (length(common) == 0) {
    warning("No matching genes found")
    return(list(scores = NULL, matched_genes = character(0)))

  } else {

    expr = t(gene_counts[match(common, rownames(gene_counts)),])
    expr = expr[,order(colnames(expr))]

    cartex_200 = cartex_200[match(common, names(cartex_200))]
    cartex_200 = cartex_200[order(names(cartex_200))]


    if(all(colnames(expr)==names(cartex_200))) {
      pca_weights = as.numeric(cartex_200)
      scores = apply(expr, 1, function(x) x %*% pca_weights)
      scores = data.frame(cbind(scores, samples=rownames(expr)))
      rownames(scores) <- rownames(expr)
      scores$cartex_score <- Zscore(scores$scores)
      scores=scores[,-1]

      if (!is.null(metadata)) {
        if (all(rownames(metadata)==rownames(expr))) {
          scores = data.frame(cbind(scores, samples=rownames(expr), metadata))
        } else {
          stop("Rownames of metadata does not match samples names of gene counts dataframe")
        }
      }

      return(scores)

    } else {
      return(paste0("genes don't match"))
    }
  }

}
