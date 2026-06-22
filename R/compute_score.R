#' Calculate PCA based phenotype specific score
#'
#' This function computes the score based on gene variance and weights
#' using principal component analysis (PCA) of a gene expression matrix
#'
#' @param gene_counts A normalized gene counts data frame with genes as rows and samples as columns.
#' @param genelist An optional character vector of gene names, Default is NULL
#' @param num_mv_genes An optional numerical variable specifying the number of most variable genes selection , default all genes from genelist or gene_counts if genelist if null, are used to compute score
#' @param pca_number An integer specifying which Principal Component to use for weights. Default is 1.
#' @return A named vector of scores.
#' @export

compute_score <- function(gene_counts, genelist=NULL, num_mv_genes=NULL,  pca_number=1) {

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

  if (!is.numeric(pca_number) || pca_number < 1 || pca_number > ncol(gene_counts)) {
    stop(paste("pc_number must be an integer between 1 and", ncol(gene_counts)))
  }

  if (is.null(num_mv_genes)) {
    if (!is.null(genelist)) {
      num_mv_genes <- length(genelist)
    } else {
      num_mv_genes <- length(gene_counts)
    }
  }

  Zscore <- function(s) {
    s <- as.numeric(s)
    z <- (s - mean(s)) / sd(s)
    return(z)
  }

  var_genes <- apply(gene_counts, 1, var)

  topXmostvariable_genes <- names(sort(var_genes, decreasing = TRUE))[1:num_mv_genes]

  pca_weights <- get_weights(gene_counts, genelist=genelist, pca_number)

  select_genecount_df <- gene_counts[match(topXmostvariable_genes, rownames(gene_counts)), ]

  select_weights_df <- pca_weights[match(topXmostvariable_genes, pca_weights$gene),]

  expr=t(select_genecount_df)

  if(all(colnames(expr)==rownames(select_weights_df))) {
    pca_weights = select_weights_df$weights
    scores = apply(expr, 1, function(x) x %*% pca_weights)
    scores = data.frame(cbind(scores, samples=rownames(expr)))
    rownames(scores) <- rownames(expr)
    scores$pca_score <- Zscore(scores$scores)
    scores=scores[,-1]
    return(scores)
  } else {
    return(paste0("genes don't match"))
  }
}
