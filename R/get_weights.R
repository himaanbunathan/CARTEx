#' Calculate Gene Weights from PCA
#'
#' This function computes weights using principal component analysis (PCA)
#' on the input gene expression matrix
#'
#' @param gene_counts A normalized gene counts data frame with genes as rows and samples as columns.
#' @param genelist An optional character vector of gene names to subset gene counts, Default is NULL
#' @param pca_number An integer specifying which Principal Component to use for weights. Default is 1.
#' @return gene weights based on PCA
#' @import stats
#' @export
get_weights <- function(gene_counts, genelist=NULL, pca_number = 1) {

  if (!is.data.frame(gene_counts)) {
    stop("Input gene_counts must be a dataframe.")
  }

  if (is.null(rownames(gene_counts))) {
    stop("The rownames of Input dataframe should be genes.")
  }

  if (!is.null(genelist)) {
    gene_counts <- gene_counts[rownames(gene_counts) %in% genelist, ]
    if (nrow(gene_counts) == 0) {
      stop("No genes in the input dataframe match the provided genelist.")
    }
  }

  # PCA and extract weights
  x <- t(gene_counts)

  pca <- prcomp(x, scale = TRUE)

  PC1_weights <- pca$rotation[, pca_number]

  df <- data.frame(gene = rownames(gene_counts), weights = PC1_weights, stringsAsFactors = FALSE)
  weights_df <- df[order(abs(df$weights), decreasing = TRUE),]

  return(weights_df)
}
