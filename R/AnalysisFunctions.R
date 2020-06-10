#' Creates clustering of the cells
#'
#' Results of this are stored as a new variabe in the object's metaData
#'
#' @param object the VISION object for which to cluster the cells
#' @return the VISION object modifed as described above
clusterCells <- function(object, tree=FALSE) {

    message("Clustering cells...", appendLF = FALSE)

    res <- object@LatentSpace

    K <- min(object@params$numNeighbors, 30)
    
    if (!tree) {
      kn <- find_knn_parallel(res, K)
    } else {
      message("Using Tree to compute clusters...\n")
      kn <- find_knn_parallel_tree(object@Tree, K)
    }
    

    cl <- louvainCluster(kn, res)

    names(cl) <- paste('Cluster', seq(length(cl)))

    # cl is list of character vector
    if (!tree) {
      cluster_variable <- "VISION_Clusters"
    } else {
      cluster_variable <- "VISION_Clusters (Tree)"
    }
    
    metaData <- object@metaData

    metaData[cluster_variable] <- factor(levels = names(cl))

    for (cluster in names(cl)) {
        metaData[cl[[cluster]], cluster_variable] <- cluster
    }

    object@metaData <- metaData

    message("completed\n")

    return(object)

}


#' create micro-clusters that reduce noise and complexity while maintaining
#' the overall signal in the data
#' @param object the VISION object for which to cluster the cells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @return the VISION with pooled cells
poolCells <- function(object, cellsPerPartition = NULL) {

    if (!is.null(cellsPerPartition)){
        object@params$micropooling$cellsPerPartition <- cellsPerPartition
    }

    if (length(object@Pools) == 0) {
        message(paste(
          "Performing micro-pooling on",
          ncol(object@exprData),
          "cells with a target pool size of",
          object@params$micropooling$cellsPerPartition
        ))
    } else {
      message("Performing micro-pooling using precomputed pools")
    }

    preserve_clusters <- NULL

    if (is.null(object@params$latentSpace[["projectionGenes"]])){
        filterInput <- object@params$latentSpace$projectionGenesMethod
    } else {
        filterInput <- object@params$latentSpace$projectionGenes
    }

    if (length(object@Pools) == 0) {
        pools <- applyMicroClustering(
            object@exprData,
            cellsPerPartition = object@params$micropooling$cellsPerPartition,
            filterInput = filterInput,
            filterThreshold = object@params$latentSpace$threshold,
            latentSpace = object@LatentSpace,
            K = object@params$numNeighbors)

        object@Pools <- pools
    }

    message("    Aggregating data using assigned pools...", appendLF = FALSE)
    pooled_cells <- poolMatrixCols(object@exprData, object@Pools)
    object@exprData <- pooled_cells

    if (hasUnnormalizedData(object)) {
        pooled_unnorm <- poolMatrixCols(object@unnormalizedData, object@Pools)
        object@unnormalizedData <- pooled_unnorm
    }

    if (!all(dim(object@LatentSpace) == c(1, 1))) {
        pooled_latent <- poolMatrixRows(object@LatentSpace, object@Pools)
        object@LatentSpace <- pooled_latent
    }

    if (hasProteinData(object)) {
        pooled_fbc <- poolMatrixRows(object@proteinData, object@Pools)
        object@proteinData <- pooled_fbc
    }

    poolMeta <- poolMetaData(object@metaData, object@Pools)
    object@metaData <- poolMeta

    if (length(object@Projections) > 0){

        newProjections <- lapply(
            object@Projections,
            function(proj) {
                new_coords <- poolMatrixRows(proj, object@Pools)
                return(new_coords)
            })

        names(newProjections) <- names(object@Projections)

        object@Projections <- newProjections

    }

    message("completed\n")

    return(object)
}

#' filter data accourding to the provided filters
#' @param object the VISION object
#' @param threshold Threshold to apply when using the 'threshold' or 'fano' projection genes filter.
#' If greater than 1, this specifies the number of cells in which a gene must be detected
#' for it to be used when computing PCA. If less than 1, this instead specifies the proportion of cells needed
#' @param num_mad Number of median absolute deviations to use when selecting highly-variable
#' genes in each mean-sorted bin of genes
#' @param projection_genes_method a method to select genes either 'Threshold' or 'Fano'
#' @param projection_genes a list of genes to use specifically.
#' If this is supplied then `projection_genes_method` is ignored
#' @return the VISION object, populated with filtered data
computeProjectionGenes <- function(object,
                       threshold = NULL,
                       num_mad = NULL,
                       projection_genes_method = NULL,
                       projection_genes = NULL) {

    if (!is.null(threshold)){
        if (threshold < 1) {
            num_samples <- ncol(object@exprData)
            threshold <- round(threshold * num_samples)
        }
        object@params$latentSpace$threshold <- threshold
    }

    if (!is.null(num_mad)){
        object@params$latentSpace$num_mad <- num_mad
    } else {
        object@params$latentSpace$num_mad <- 2
    }

    if (!is.null(projection_genes)){
        object@params$latentSpace$projectionGenes <- projection_genes
    } else {
        object@params$latentSpace$projectionGenes <- NULL
    }

    if (!is.null(projection_genes_method)){
        object@params$latentSpace$projectionGenesMethod <- projection_genes_method
    }

    message("Determining projection genes...")

    if (is.null(object@params$latentSpace[["projectionGenes"]])){

        exprData <- matLog2(object@exprData)
        projection_genes <- applyFilters(
                    exprData,
                    object@params$latentSpace$projectionGenesMethod,
                    object@params$latentSpace$threshold,
                    object@params$latentSpace$num_mad)

        if (length(projection_genes) == 0){
            stop(
                sprintf("Filtering with (projection_genes=\"%s\", threshold=%i) results in 0 genes\n  Set a lower threshold and re-run",
                    object@params$latentSpace$projectionGenesMethod, object@params$latentSpace$threshold)
                )
        }
    } else {
        projection_genes <- intersect(
            object@params$latentSpace$projectionGenes, rownames(object@exprData))

        if (length(projection_genes) == 0){
            stop("Supplied list of genes in `projection_genes` does not match any rows of expression data")
        } else {
            message(
                sprintf("    Using supplied list of genes: Found %i/%i matches",
                    length(projection_genes), length(object@params$latentSpace$projectionGenes)
                    )
                )
        }
    }

    message()

    object@params$latentSpace$projectionGenes <- projection_genes


    return(object)
}


#' Add signatures to VISION object
#'
#'
#' @param object the VISION object
#' @param signatures list of file paths to signature files (.gmt or .txt) or
#' Signature objects.  See the createGeneSignature(...) method for information
#' on creating Signature objects.
#' @param min_signature_genes Signature that match less than this number of genes in the
#' supplied expression matrix are removed.
#' @param sig_gene_threshold Proportion of cells that a gene must be detected in (nonzero)
#' to be used in signature score calculations.
#' @return the VISION object, with the @sigData slot updated
#' @export
addSignatures <- function(object, signatures, min_signature_genes=5, sig_gene_threshold=.01) {

    if (is.list(signatures)) {
        sigs <- lapply(signatures, function(sig){
                   if (is(sig, "Signature")){
                       return(sig)
                   } else {
                       return(readSignaturesInput(sig))
                   }
        })

        if (length(sigs) > 0){
            sigs <- do.call(c, sigs)
        }

        names(sigs) <- vapply(sigs, function(x){x@name}, "")

    } else if (is.character(signatures)) {
        sigs <- readSignaturesInput(signatures)
    } else {
        stop("signatures must be paths to signature files or list of
            Signature objects")
    }

    sigs <- processSignatures(sigs, object@exprData, min_signature_genes, sig_gene_threshold)

    object@sigData <- c(object@sigData, sigs)

    return(object)
}


#' calculate signature scores
#'
#' For each signature-cell pair, compute a score that captures the level of
#' correspondence between the cell and the signature.
#'
#' @param object the VISION object
#' @param sig_norm_method Method to apply to normalize the expression matrix
#' before calculating signature scores. Valid options are:
#' "znorm_columns" (default), "none", "znorm_rows", "znorm_rows_then_columns",
#' or "rank_norm_columns"
#' @param sig_gene_importance whether or not to rank each gene's contribution to
#' the overall signature score.  Default = TRUE.  This is used for inspecting
#' genes in a signature in the output report
#' @return the VISION object, with the @SigScores and @SigGeneImportance slots populated
#' @export
calcSignatureScores <- function(
    object, sig_norm_method = NULL, sig_gene_importance = TRUE) {

    message("Evaluating signature scores on cells...\n")

    ## override object parameters
    if (!is.null(sig_norm_method)) object@params$signatures$sigNormMethod <- sig_norm_method

    if (length(object@sigData) == 0) {
        sigScores <- matrix(nrow = ncol(object@exprData), ncol = 0,
                            dimnames = list(colnames(object@exprData), NULL)
                            )
        object@SigScores <- sigScores
        object@SigGeneImportance <- list()
        return(object)
    }

    normExpr <- getNormalizedCopySparse(
        object@exprData,
        object@params$signatures$sigNormMethod
    )

    sigScores <- batchSigEvalNorm(object@sigData, normExpr)

    if (sig_gene_importance) {

        if (is(object@exprData, "sparseMatrix")) {
            sigGeneImportance <- evalSigGeneImportanceSparse(
                sigScores, object@sigData, normExpr
            )
        } else {
            normExprDense <- getNormalizedCopy(
                object@exprData,
                object@params$signatures$sigNormMethod
            )
            sigGeneImportance <- evalSigGeneImportance(
                sigScores, object@sigData, normExprDense
            )
        }

    } else {
        sigGeneImportance <- list()
    }

    object@SigScores <- sigScores
    object@SigGeneImportance <- sigGeneImportance

    return(object)
}


#' Calculate gene-signature importance
#'
#' For each signature, the contribution of each gene to the signature score
#' is evaluated by calculating the covariance between signature scores and expression
#' The correlation of genes with a negative sign in the signature are inverted.
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colSds
#' @importFrom matrixStats rowSds
#' @importFrom stats setNames
#'
#' @param sigScores matrix of signature scores (Cells x Signatures)
#' @param sigData list of Signature objects
#' @param normExpr output from calling getNormalizedCopySparse
#' @return named list with one item per signature.  Values are named vectors
#' of each gene's association (covariance) with the signature.
evalSigGeneImportance <- function(sigScores, sigData, normExpr){

    message("Evaluating signature-gene importance...\n")

    if (length(sigData) == 0) {
        return(list())
    }

    if (length(sigScores) <= 1){
        stop(
            sprintf("Signature scores have not yet been computed.  `calcSignatureScores` must be run before running `evalSigGeneImportance`")
            )
    }

    # Center each column of sigScores first

    mu <- colMeans(sigScores)

    sigScores <- t(sigScores)
    sigScores <- (sigScores - mu)
    sigScores <- t(sigScores)

    # Center each row of normExpr
    mu <- rowMeans(normExpr)

    normExpr <- (normExpr - mu)

    # Compute Covariances
    sigGene <- function(signame) {
        sigdata <- sigData[[signame]]

        genes <- sigdata@sigDict

        sigvals <- sigScores[, signame]

        geneIndices <- match(names(genes), rownames(normExpr))

        corr <- sigGeneInner(sigvals, normExpr, geneIndices)

        names(corr) <- names(genes)

        corr <- corr * genes

        return(corr)
    }

    sigs <- colnames(sigScores)
    res <- pbmclapply(setNames(sigs, sigs), sigGene)

    return(res)
}


#' Calculate Gene-Signature Importance
#'
#' For each signature, the contribution of each gene to the signature score
#' is evaluated by calculating the covariance between signature scores and expression
#' The correlation of genes with a negative sign in the signature are inverted.
#'
#' This version is made to avoid inflating sparse matrices
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colSds
#' @importFrom matrixStats rowSds
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#' @importFrom stats setNames
#'
#' @param object the VISION object
#' @return the VISION object, with SigGeneImportance slot populated
#' @param sigScores matrix of signature scores (Cells x Signatures)
#' @param sigData list of Signature objects
#' @param normExpr output from calling getNormalizedCopySparse
#' @return named list with one item per signature.  Values are named vectors
#' of each gene's association (covariance) with the signature.
evalSigGeneImportanceSparse <- function(sigScores, sigData, normExpr){

    message("Evaluating signature-gene importance...\n")

    if (length(sigData) == 0) {
        return(list())
    }

    if (length(sigScores) <= 1){
        stop(
            sprintf("Signature scores have not yet been computed.  `calcSignatureScores` must be run before running `evalSigGeneImportance`")
            )
    }

    # Center each column of sigScores first
    mu <- colMeans(sigScores)

    sigScores <- t(sigScores)
    sigScores <- (sigScores - mu)
    sigScores <- t(sigScores)

    # Precompute some matrices we'll need later
    NGenes <- nrow(normExpr@data)
    NCells <- ncol(normExpr@data)
    Cog <- Matrix(1, ncol = 1, nrow = NGenes)
    Coc <- Matrix(normExpr@colOffsets, nrow = 1)
    Cs <- Diagonal(x = normExpr@colScaleFactors)
    Roc <- Matrix(1, nrow = 1, ncol = NCells)
    C1 <- t(Roc)

    RM <- normExpr@data %*% (Cs %*% C1) + Cog %*% (Coc %*% (Cs %*% C1))
    RM <- RM / NCells
    RM <- RM[, 1]

    # Compute Covariances
    sigGene <- function(signame) {
        sigdata <- sigData[[signame]]

        genes <- sigdata@sigDict

        S <- sigScores[, signame, drop = F]
        geneIndices <- rownames(normExpr@data) %in% names(genes)

        E <- normExpr@data[geneIndices, , drop = FALSE]

        Rog <- Matrix(RM[geneIndices], ncol = 1)
        Cog <- Matrix(1, ncol = 1, nrow = length(genes))

        geneCov <- E %*% Cs %*% S - Rog %*% (Roc %*% S) + Cog %*% (Coc %*% (Cs %*% S))
        geneCov <- geneCov / (NCells - 1)
        geneCov <- geneCov[, 1]
        geneCov <- geneCov[names(genes)]

        geneCov <- geneCov * genes # invert sign for negative genes

        return(geneCov)
    }

    sigs <- colnames(sigScores)
    res <- pbmclapply(setNames(sigs, sigs), sigGene)

    return(res)
}


#' Computes the latent space of the expression matrix using PCA
#'
#' @param object the VISION object for which compute the latent space
#' @param projection_genes character vector of gene names to use for PCA
#' @param projection_genes_method name of filtering method. Either 'threshold' or 'fano'(default)
#' @param filterThreshold Threshold to apply when using the 'threshold' or 'fano' projection genes filter.
#' If greater than 1, this specifies the number of cells in which a gene must be detected
#' for it to be used when computing PCA. If less than 1, this instead specifies the proportion of cells needed
#' @param filterNumMad Number of median absolute deviations to use when selecting highly-variable
#' genes in each mean-sorted bin of genes
#' @param num_PCs the number of principal components to retain
#' @param perm_wPCA If TRUE, apply permutation wPCA to determine significant
#' number of components. Default is FALSE.
#' @return the VISION with @latentSpace slot populated
#' @export
computeLatentSpace <- function(
    object, projection_genes = NULL,
    filterThreshold = .05, filterNumMad = 2,
    projection_genes_method = NULL,
    num_PCs = 30, perm_wPCA = NULL) {

    message("Computing a latent space for expression data...\n")

    if (!is.null(projection_genes)) {
        object@params$latentSpace$projectionGenes <- projection_genes
    }
    if (!is.null(projection_genes_method)) {
        object@params$latentSpace$projectionGenesMethod <- projection_genes_method
    }

    if (is.null(object@params$latentSpace[["projectionGenes"]])) {
        object <- computeProjectionGenes(
            object,
            threshold = filterThreshold,
            num_mad = filterNumMad,
            projection_genes_method =
                object@params$latentSpace$projectionGenesMethod
        )
    } else {
        object <- computeProjectionGenes(
            object,
            projection_genes = object@params$latentSpace$projectionGenes
        )
    }

    if (!is.null(perm_wPCA)) object@params$latentSpace$permPCA <- perm_wPCA
    object@params$latentSpace$numPCs <- num_PCs

    expr <- object@exprData
    projection_genes <- object@params$latentSpace[["projectionGenes"]]
    perm_wPCA <- object@params$latentSpace$permPCA

    if (!is.null(projection_genes)) {
        exprData <- expr[projection_genes, , drop = FALSE]
    } else {
        exprData <- expr
    }

    exprData <- matLog2(exprData)

    if (perm_wPCA) {
        res <- applyPermutationWPCA(exprData, components = num_PCs)
        pca_res <- res[[1]]
    } else {
        res <- applyPCA(exprData, maxComponents = num_PCs)
        pca_res <- res[[1]]
    }

    object@LatentSpace <- t(pca_res)
    colnames(object@LatentSpace) <- paste0("PC ", seq_len(ncol(object@LatentSpace)))
    object@params$latentSpace$name <- "PCA"
    return(object)
}


#' Add a latent space computed using an external method
#'
#' @param object the VISION object for which compute the latent space
#' @param coordinates matrix with latent space coordinates (cells x dimensions)
#' @param name a label for the latent space (e.g. "PCA" or "scVI")
#' @return the VISION with @latentSpace slot populated
addLatentSpace <- function(object, coordinates, name) {

    if (is.data.frame(coordinates)){
        coordinates <- data.matrix(coordinates)
    }

    if (is.null(rownames(coordinates))) {
        if (nrow(coordinates) != ncol(object@exprData)) {
            stop("Supplied coordinates must of number of rows equal to number of cells in expression matrix")
        }

        rownames(coordinates) <- colnames(object@exprData)
    } else {
        sample_names <- colnames(object@exprData)
        common <- intersect(sample_names, rownames(coordinates))

        if (length(common) != nrow(coordinates)){
            stop("Supplied coordinates for coordinates must have rowlabels that match sample/cell names")
        }
        coordinates <- coordinates[colnames(object@exprData), , drop = FALSE]
    }

    if (is.null(colnames(coordinates))) {
        colnames(coordinates) <- paste0(name, "-", seq_len(ncol(coordinates)))
    }

    object@LatentSpace <- coordinates
    object@params$latentSpace$name <- name
    return(object)
}


#' generate projections
#'
#' Generates 2-dimensional representations of the expression matrix
#' Populates the 'Projections' slot on the VISION object
#'
#' @importFrom stats setNames
#'
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
generateProjections <- function(object) {
  message("Projecting data into 2 dimensions...")

  # Some projection methods operate on the full expression matrix
  # If using one of these, we need to compute 'projection_genes'
  projection_methods <- object@params$projectionMethods
  if ("ICA" %in% projection_methods || "RBFPCA" %in% projection_methods) {
      object <- computeProjectionGenes(object)
  }

  projections <- generateProjectionsInner(object@exprData,
                                     object@LatentSpace,
                                     projection_genes = object@params$latentSpace[["projectionGenes"]],
                                     projection_methods = object@params$projectionMethods,
                                     K = object@params$numNeighbors)

  # Add already-input projections
  for (proj in names(object@Projections)){
      projections[[proj]] <- object@Projections[[proj]]
  }

  # Make sure all projections have column names
  n <- names(projections)
  projections <- lapply(setNames(n, n), function(pname){
      proj <- projections[[pname]]
      if (is.null(colnames(proj))){
          colnames(proj) <- paste0(pname, "-", seq_len(ncol(proj)))
      }
      return(proj)
  })

  object@Projections <- projections

  message("")

  return(object)
}


#' Adds a UMAP projection
#'
#' @param object the VISION object
#' @param K Number of neighbors to use in UMAP projection.
#' @param name label to use for this projection
#' @param source coordinates to use to compute tSNE
#'   Choices are either "LatentSpace" (default) or "Proteins"
#'
#' @return VISION object with the projection added to @Projections slot
#' @export
addUMAP <- function(object, K = object@params$numNeighbors,
                    name = "UMAP", source = "LatentSpace") {

    if (!requireNamespace("uwot", quietly = TRUE)){
        stop("Package \"uwot\" needed to run UMAP.  Please install it using:\n\n   devtools::install_github(\"jlmelville/uwot\")\n\n",
            call. = FALSE)
    }

    if (source == "LatentSpace") {
        data <- object@LatentSpace
    } else if (source == "Proteins") {
        data <- object@proteinData
    } else {
        stop("Invalid 'source' parameter.  Can be either 'LatentSpace' or 'Proteins'")
    }

    n_workers <- getOption("mc.cores")
    n_workers <- if (is.null(n_workers)) 2 else n_workers
	res <- uwot::umap(
        data, n_neighbors = K,
        n_threads = n_workers, ret_nn = T
    )
	res <- res$embedding

	rownames(res) <- rownames(data)
    colnames(res) <- paste0(name, "-", seq_len(ncol(res)))
    object@Projections[[name]] <- res

    return(object)
}


#' Adds a tSNE projection
#'
#' @importFrom Rtsne Rtsne
#'
#' @param object the VISION object
#' @param perplexity parameter for tSNE
#' @param name label to use for this projection
#' @param source coordinates to use to compute tSNE
#'   Choices are either "LatentSpace" (default) or "Proteins"
#' @return VISION object with the projection added to @Projections slot
#' @export
addTSNE <- function(object, perplexity = 30, name = "tSNE", source = "LatentSpace") {

    if (source == "LatentSpace") {
        data <- object@LatentSpace
    } else if (source == "Proteins") {
        data <- object@proteinData
    } else {
        stop("Invalid 'source' parameter.  Can be either 'LatentSpace' or 'Proteins'")
    }

    res <- Rtsne(
        data, dims = 2, max_iter = 800, perplexity = perplexity,
        check_duplicates = FALSE, pca = FALSE)
    res <- res$Y
    rownames(res) <- rownames(data)
    colnames(res) <- paste0(name, "-", seq_len(ncol(res)))

    object@Projections[[name]] <- res

    return(object)

}

#' Compute local correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
#' @export
analyzeLocalCorrelations <- function(object, tree=FALSE) {

  signatureBackground <- generatePermutationNull(
      object@exprData, object@sigData, num = 3000
  )

  normExpr <- getNormalizedCopySparse(
      object@exprData,
      object@params$signatures$sigNormMethod)

  message("Computing KNN Cell Graph in the Latent Space...\n")
  if (!tree) {
    weights <- computeKNNWeights(object@LatentSpace, object@params$numNeighbors)
  } else {
    message("Using Tree to compute neighbors...\n")
    weights <- computeKNNWeights(object@Tree, object@params$numNeighbors)
  }
 
  message("Evaluating local consistency of signatures in latent space...\n")

  sigConsistencyScores <- sigConsistencyScores(
                                weights,
                                object@SigScores,
                                object@metaData,
                                signatureBackground,
                                normExpr)

  if (hasProteinData(object)) {
      fbcs <- fbConsistencyScores(weights, object@proteinData)
  } else {
      fbcs <- NULL
  }

  message("Clustering signatures...\n")
  sigClusters <- clusterSignatures(object@SigScores,
                                   object@metaData,
                                   sigConsistencyScores,
                                   clusterMeta = object@params$micropooling$pool)

  metaConsistencyScores <- sigConsistencyScores[
      colnames(object@metaData), , drop = FALSE
  ]

  sigConsistencyScores <- sigConsistencyScores[
      colnames(object@SigScores), , drop = FALSE
  ]


  LocalAutocorrelation <- list(
      "Signatures" = sigConsistencyScores,
      "Meta" = metaConsistencyScores,
      "Proteins" = fbcs,
      "Clusters" = sigClusters
  )

  object@LocalAutocorrelation <- LocalAutocorrelation

  return(object)
}


#' Compute trajectory correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
analyzeTrajectoryCorrelations <- function(object) {

  signatureBackground <- generatePermutationNull(
      object@exprData, object@sigData, num = 3000
  )

  normExpr <- getNormalizedCopySparse(
      object@exprData,
      object@params$signatures$sigNormMethod)

  message("Computing KNN Cell Graph in the Trajectory Model...\n")

  weights <- computeKNNWeights(object@LatentTrajectory, object@params$numNeighbors)

  message("Evaluating local consistency of signatures within trajectory model...\n")

  sigVTreeProj <- sigConsistencyScores(weights,
                                       object@SigScores,
                                       object@metaData,
                                       signatureBackground,
                                       normExpr)

  if (hasProteinData(object)) {
      fbcs <- fbConsistencyScores(weights, object@proteinData)
  } else {
      fbcs <- NULL
  }

  message("Clustering signatures...\n")
  sigTreeClusters <- clusterSignatures(object@SigScores,
                                       object@metaData,
                                       sigVTreeProj,
                                       clusterMeta = object@params$micropooling$pool)

  metaConsistencyScores <- sigVTreeProj[
      colnames(object@metaData), , drop = FALSE
  ]

  sigConsistencyScores <- sigVTreeProj[
      colnames(object@SigScores), , drop = FALSE
  ]

  TrajectoryAutocorrelation <- list(
      "Signatures" = sigConsistencyScores,
      "Meta" = metaConsistencyScores,
      "Proteins" = fbcs,
      "Clusters" = sigTreeClusters
  )

  object@TrajectoryAutocorrelation <- TrajectoryAutocorrelation

  return(object)
}


#' Compute Ranksums Test, for all factor meta data.  One level vs all others
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colRanks
#' @importFrom stats setNames
#' @param object the VISION object
#' @param variables which columns of the meta-data to use for comparisons
#' @return the VISION object with the @ClusterComparisons slot populated
#' @export
clusterSigScores <- function(object, variables = "All") {

    message("Computing differential signature tests...\n")

    sigScores <- object@SigScores
    metaData <- object@metaData

    metaData <- metaData[rownames(sigScores), , drop = FALSE]

    if (variables == "All") {
        # Determine which metaData we can run on
        # Must be a factor with at least 20 levels
        clusterMeta <- vapply(colnames(metaData), function(x) {
                scores <- metaData[[x]]
                if (!is.factor(scores)){
                    return("")
                }
                if (length(levels(scores)) > 50){
                    return("")
                }
                if (length(unique(scores)) == 1){
                    return("")
                }
                return(x)
            }, FUN.VALUE = "")
        clusterMeta <- clusterMeta[clusterMeta != ""]
    } else {
        if (!all(variables %in% colnames(metaData))) {
            stop("Supplied variable names must be column names of object@metaData")
        }
        clusterMeta <- setNames(variables, variables)
    }

    ClusterComparisons <- list()

    # Comparisons for Signatures
    if (ncol(sigScores) > 0){
        sigScoreRanks <- colRanks(sigScores,
                                  preserveShape = TRUE,
                                  ties.method = "average")
        dimnames(sigScoreRanks) <- dimnames(sigScores)
    } else {
        sigScoreRanks <- sigScores
    }

    out <- pbmclapply(clusterMeta, function(variable){
        values <- metaData[[variable]]
        var_levels <- levels(values)

        result <- lapply(var_levels, function(var_level){
            cluster_ii <- which(values == var_level)

            r1 <- matrix_wilcox(sigScoreRanks, cluster_ii,
                                check_na = FALSE, check_ties = FALSE)

            pval <- r1$pval
            stat <- r1$stat
            fdr <- p.adjust(pval, method = "BH")
            out <- data.frame(
                stat = stat, pValue = pval, FDR = fdr
            )
            return(out)
        })

        names(result) <- var_levels
        result <- result[order(var_levels)]

        return(result)
    }, mc.cores = 1)

    ClusterComparisons[["Signatures"]] <- out

    # Comparisons for Meta data
    # Split meta into numeric and factor
    numericMeta <- vapply(seq_len(ncol(metaData)),
                          function(i) is.numeric(metaData[[i]]),
                          FUN.VALUE = TRUE)
    numericMeta <- metaData[, numericMeta, drop = F]
    numericMeta <- as.matrix(numericMeta)

    factorMeta <- vapply(seq_len(ncol(metaData)),
                          function(i) is.factor(metaData[[i]]),
                          FUN.VALUE = TRUE)
    factorMeta <- metaData[, factorMeta, drop = F]

    if (ncol(numericMeta) > 0){
        numericMetaRanks <- colRanks(numericMeta,
                                     preserveShape = TRUE,
                                     ties.method = "average")
        dimnames(numericMetaRanks) <- dimnames(numericMeta)
    } else {
        numericMetaRanks <- numericMeta
    }

    out <- pbmclapply(clusterMeta, function(variable){
        values <- metaData[[variable]]
        var_levels <- levels(values)

        result <- lapply(var_levels, function(var_level){
            cluster_ii <- which(values == var_level)

            r2 <- matrix_wilcox(numericMetaRanks, cluster_ii,
                                check_na = TRUE, check_ties = TRUE)

            r3 <- matrix_chisq(factorMeta, cluster_ii)

            pval <- c(r2$pval, r3$pval)
            stat <- c(r2$stat, r3$stat)
            fdr <- p.adjust(pval, method = "BH")
            out <- data.frame(
                stat = stat, pValue = pval, FDR = fdr
            )
            return(out)
        })

        names(result) <- var_levels
        result <- result[order(var_levels)]

        return(result)
    }, mc.cores = 1)

    ClusterComparisons[["Meta"]] <- out

    # Comparisons for Proteins
    if (hasProteinData(object)){
        fbData <- t(object@proteinData)

        # gather jobs
        jobs <- lapply(unname(clusterMeta), function(variable) {
            lapply(levels(metaData[[variable]]), function(level) {
                return(c(variable, level))
            })
        })
        jobs <- do.call(c, jobs)

        results <- pbmclapply(jobs, function(job){
            variable <- job[1]
            var_level <- job[2]
            values <- metaData[[variable]]
            cluster_ii <- which(values == var_level)
            not_cluster_ii <- which(values != var_level)

            rr <- matrix_wilcox_cpp(
                fbData, cluster_ii, not_cluster_ii, jobs = 1)

            pval <- rr$pval
            stat <- rr$AUC
            fdr <- p.adjust(pval, method = "BH")
            out <- data.frame(
                stat = stat, pValue = pval, FDR = fdr,
                row.names = rownames(rr)
            )
            return(list(job, out))
        })

        out_fb <- list()
        for (res in results){
            variable <- res[[1]][1]
            var_level <- res[[1]][2]
            df <- res[[2]]
            if (is.null(out_fb[[variable]])){
                out_fb[[variable]] <- list()
            }
            out_fb[[variable]][[var_level]] <- df
        }

        ClusterComparisons[["Proteins"]] <- out_fb
    } else {
        ClusterComparisons[["Proteins"]] <- NULL
    }

    object@ClusterComparisons <- ClusterComparisons

    return(object)

}

#' Compute pearson correlation between signature scores and components of the Latent Space
#'
#' Enables the LCAnnotator portion of the output report
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats cor.test
#' @param object the VISION object
#' @return object with the @LCAnnotatorData slot populated
#' @export
annotateLatentComponents <- function(object){

  message("Computing correlations between signatures and latent space components...\n")
  sigMatrix <- object@SigScores
  metaData <- object@metaData
  latentSpace <- object@LatentSpace

  ## combined gene signs and numeric meta variables

  numericMetaVars <- vapply(colnames(metaData),
                            function(x) is.numeric(metaData[[x]]),
                            FUN.VALUE = TRUE)

  numericMeta <- metaData[, numericMetaVars, drop = FALSE]
  numericMeta <- numericMeta[rownames(sigMatrix), , drop = FALSE]


  computedSigMatrix <- cbind(sigMatrix, numericMeta)

  pearsonCorr <- pbmclapply(seq_len(ncol(computedSigMatrix)), function(i) {
      ss <- computedSigMatrix[, i];

      ls_col_cor <- apply(latentSpace, 2, function(pc){
           suppressWarnings({
               pc_result <- cor.test(ss, pc)
           })
           if (is.na(pc_result$estimate)) {  # happens if std dev is 0 for a sig
               return(0)
           } else {
               return(pc_result$estimate)
           }
      })

      return(ls_col_cor)
  })

  if (length(pearsonCorr) > 0) {
      pearsonCorr <- do.call(rbind, pearsonCorr)
  } else {
      pearsonCorr <- matrix(
          nrow = ncol(computedSigMatrix),
          ncol = ncol(latentSpace)
      )
  }
  rownames(pearsonCorr) <- colnames(computedSigMatrix)
  colnames(pearsonCorr) <- colnames(latentSpace)

  if (hasProteinData(object)) {
      proteinData <- object@proteinData
      pearsonCorrProteins <- pbmclapply(
          seq_len(ncol(proteinData)), function(i) {
          ss <- proteinData[, i];

          ls_col_cor <- apply(latentSpace, 2, function(pc){
               suppressWarnings({
                   pc_result <- cor.test(ss, pc)
               })
               if (is.na(pc_result$estimate)) {  # happens if std dev is 0 for a protein
                   return(0)
               } else {
                   return(pc_result$estimate)
               }
          })

          return(ls_col_cor)
      })

      pearsonCorrProteins <- do.call(rbind, pearsonCorrProteins)
      rownames(pearsonCorrProteins) <- colnames(proteinData)
      colnames(pearsonCorrProteins) <- colnames(latentSpace)
  } else {
      pearsonCorrProteins <- NULL
  }


  pcaAnnotData <- LCAnnotatorData(
      pearsonCorr = pearsonCorr, pearsonCorrProteins = pearsonCorrProteins
  )
  object@LCAnnotatorData <- pcaAnnotData

  return(object)
}
