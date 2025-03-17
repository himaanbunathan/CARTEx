#' Apply CARTEx Score to Each Cell in a Seurat Object
#'
#' This function computes the CARTEx score for each cell in a Seurat object 
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param assay The assay in the Seurat object to use (default = "RNA").
#' @return A Seurat object with CARTEx scores added to meta data.
#' @import Seurat
#' @export
#'
single_cell_cartex <- function(seurat_obj, assay = "RNA") {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat v5 object.")
  }
  
  Zscore <- function(s) {
    s <- as.numeric(s)
    z <- (s - mean(s)) / sd(s)
    return(z)
  }
  
  if (!is.null(seurat_obj[[assay]]@scale.data) && nrow(seurat_obj[[assay]]@scale.data) > 0) {
    message("Using scaled data for CARTEx computation.")
    gene_counts <- GetAssayData(seurat_obj, layer = "scale.data", assay = assay)
    
  } else if (!is.null(seurat_obj[[assay]]@data) && nrow(seurat_obj[[assay]]@data) > 0) {
    message("Using normalized data for CARTEx computation.")
    gene_counts <- GetAssayData(seurat_obj, layer = "data", assay = assay)
    
  } else {
    message("Warning: Using raw counts for CARTEx computation. This may affect results.Consider normalizing your data before running CARTEx.")
    gene_counts <- GetAssayData(seurat_obj, layer = "counts", assay = assay)
  }
  
  # Check if gene names follow Ensembl ID format (start with "ENSG")
  if (all(grepl("^ENSG", rownames(gene_counts)))) {
    message("Detected Ensembl gene IDs. Loading Ensembl-based weight file...")
    precomputed_weights <- system.file("extdata", "weights_ensemble.rds", package = "CARTEx")
    if (precomputed_weights == "") {
      stop("precomputed weights file not found.")
    }
    
  } else {
    message("Detected gene symbols. Loading gene-symbol-based weight file...")
    precomputed_weights <- system.file("extdata", "weights.rds", package = "CARTEx")
    if (precomputed_weights == "") {
      stop("precomputed weights file not found.")
    }
  }
  
  cartex_200 <- readRDS(precomputed_weights)
  
  common_genes <- intersect(names(cartex_200), rownames(gene_counts))
  message("Number of matched genes: ", length(common_genes))
  
  if (length(common_genes) == 0) {
    warning("No matching genes found between CARTEx gene list and Seurat object features.")
    
  } else {
    print(paste0("Num of genes used for cartex score = ", length(common_genes)))
    
    expr <- t(data.frame(gene_counts[common_genes, ]))
    cartex_200 <-  cartex_200[common_genes] #cartex_weights
    
    if(all(colnames(expr)==names(cartex_200))) {
      pca_weights = as.numeric(cartex_200)
      scores = apply(expr, 1, function(x) x %*% pca_weights)
      scores = data.frame(cbind(scores, samples=rownames(expr)))
      rownames(scores) <- rownames(expr)
      scores$cartex_score <- Zscore(scores$scores)
      scores=scores[,-1]
    }

    seurat_cells <- colnames(seurat_obj)
    score_cells <- rownames(scores)
    
    print(paste0("Num of cell/barcodes in seurat object = ", length(colnames(seurat_obj))))
    print(paste0("Num of cell/barcodes used for cartex score = ", length(rownames(scores))))
    
    if (!all(seurat_cells %in% score_cells)) {
      warning("Mismatch between Seurat cell barcodes and CARTEx scores. Aligning scores...")
    }
    
    matched_scores <- rep(NA, length(seurat_cells))
    names(matched_scores) <- seurat_cells
    matched_scores[score_cells] <- scores$cartex_score
    
    seurat_obj$cartex_score <- matched_scores
    
  }
  return(seurat_obj)
}