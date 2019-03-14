#' Creates clustering of the cells
#'
#' Results of this are stored as a new variabe in the object's metaData
#'
#' @param object the VISION object for which to cluster the cells
#' @return the VISION object modifed as described above
clusterCells <- function(object) {

    message("Clustering cells...", appendLF = FALSE)

    res <- object@latentSpace

    n_workers <- getOption("mc.cores")
    n_workers <- if (is.null(n_workers)) 2 else n_workers

    kn <- ball_tree_knn(res,
                        min(round(sqrt(nrow(res))), 30),
                        n_workers)

    cl <- louvainCluster(kn, res)

    names(cl) <- paste('Cluster', seq(length(cl)))

    # cl is list of character vector
    cluster_variable <- "VISION_Clusters"
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
poolCells <- function(object,
                      cellsPerPartition=object@cellsPerPartition) {
    object@cellsPerPartition <- cellsPerPartition

    if (length(object@pools) == 0) {
        message(paste(
          "Performing micro-pooling on",
          ncol(object@exprData),
          "cells with a target pool size of",
          object@cellsPerPartition
        ))
    } else {
      message("Performing micro-pooling using precomputed pools")
    }

    preserve_clusters <- NULL

    if (length(object@pools) == 0) {
        pools <- applyMicroClustering(object@exprData,
                                              cellsPerPartition = object@cellsPerPartition,
                                              filterInput = object@projection_genes,
                                              filterThreshold = object@threshold,
                                              latentSpace = object@latentSpace)

        object@pools <- pools
    }

    message("    Aggregating data using assigned pools...", appendLF = FALSE)
    pooled_cells <- poolMatrixCols(object@exprData, object@pools)
    object@exprData <- pooled_cells

    if (hasUnnormalizedData(object)) {
        pooled_unnorm <- poolMatrixCols(object@unnormalizedData, object@pools)
        object@unnormalizedData <- pooled_unnorm
    }

    if (!all(dim(object@latentSpace) == c(1, 1))) {
        pooled_latent <- poolMatrixRows(object@latentSpace, object@pools)
        object@latentSpace <- pooled_latent
    }

    poolMeta <- poolMetaData(object@metaData, object@pools)
    object@metaData <- poolMeta

    if (length(object@inputProjections) > 0){

        newInputProjections <- lapply(
            object@inputProjections,
            function(proj) {
                new_coords <- poolMatrixRows(proj, object@pools)
                return(new_coords)
            })

        names(newInputProjections) <- names(object@inputProjections)

        object@inputProjections <- newInputProjections

    }

    message("completed\n")

    return(object)
}

#' filter data accourding to the provided filters
#' @param object the VISION object
#' @param threshold threshold to apply for the threshold filter
#' @param projection_genes either a list of genes or a method to select genes
#' @return the VISION object, populated with filtered data
filterData <- function(object,
                       threshold=object@threshold,
                       projection_genes=object@projection_genes) {

    object@projection_genes <- projection_genes
    object@threshold <- threshold

    message("Determining projection genes...")

    if (length(object@projection_genes) == 1){

        exprData <- matLog2(object@exprData)
        projection_genes <- applyFilters(
                    exprData,
                    object@threshold,
                    object@projection_genes)

        if (length(projection_genes) == 0){
            stop(
                sprintf("Filtering with (projection_genes=\"%s\", threshold=%i) results in 0 genes\n  Set a lower threshold and re-run",
                    object@projection_genes, object@threshold)
                )
        }
    } else {
        projection_genes <- intersect(
            object@projection_genes, rownames(object@exprData))

        if (length(projection_genes) == 0){
            stop("Supplied list of genes in `projection_genes` does not match any rows of expression data")
        } else {
            message(
                sprintf("    Using supplied list of genes: Found %i/%i matches",
                    length(projection_genes), length(object@projection_genes)
                    )
                )
        }
    }

    message()

    object@projection_genes <- projection_genes


    return(object)
}

#' Calculate weights based on the false negative rates, computed using the
#' provided housekeeping genes
#' @param object the VISION object
#' @param nomodel (optional) if TRUE, no fnr curve calculated and all weights
#' equal to 1. Else FNR and weights calculated.
#' @return the VISION object with populated weights slot
calcWeights <- function(object,
                        nomodel=object@nomodel) {
    object@nomodel <- nomodel
    clustered <- (ncol(object@exprData) > 15000 || object@pool)
    if (!clustered && !object@nomodel) {
        message("Computing weights from False Negative Rate curves...")
        falseneg_out <- createFalseNegativeMap(object@exprData,
                                               object@housekeepingData)

        object@weights <- computeWeights(falseneg_out[[1]], falseneg_out[[2]],
                                         object@exprData)

        rownames(object@weights) <- rownames(object@exprData)
        colnames(object@weights) <- colnames(object@exprData)

    }

    return(object)
}

#' calculate signature scores
#'
#' For each signature-cell pair, compute a score that captures the level of
#' correspondence between the cell and the signature.
#' To estimate significance of these scores, a set of random gene signatures is
#' generated to create a null distribution
#'
#' @param object the VISION object
#' @param sigData a list of Signature objects for which to compute the scores
#' @param metaData a list of existing cell-signatures
#' @param sig_norm_method (optional) Method to apply to normalize the expression
#' matrix before calculating signature scores
#' @param sig_score_method the scoring method to use
#' @return the VISION object, with signature score slots populated
calcSignatureScores <- function(object,
                                sigData=object@sigData,
                                metaData=object@metaData,
                                sig_norm_method=object@sig_norm_method,
                                sig_score_method=object@sig_score_method) {

    message("Evaluating signature scores on cells...\n")

    ## override object parameters
    if(!is.null(sigData)) object@sigData <- sigData
    if(!is.null(metaData)) object@metaData <- metaData
    object@sig_norm_method <- sig_norm_method
    object@sig_score_method <- sig_score_method

    if (length(object@sigData) == 0) {
        sigScores <- matrix(nrow=ncol(object@exprData), ncol=0,
                            dimnames = list(colnames(object@exprData), NULL)
                            )
        object@sigScores <- sigScores
        return(object)
    }

    normExpr <- getNormalizedCopy(object@exprData, object@sig_norm_method)

    sigScores <- batchSigEval(object@sigData, object@sig_score_method,
                              normExpr, object@weights)

    object@sigScores <- sigScores

    return(object)
}

#' Computes the latent space of the expression matrix using PCA
#'
#' @param object the VISION object for which compute the latent space
#' @param projection_genes character vector of gene names to use for projections
#' @param perm_wPCA If TRUE, apply permutation wPCA to determine significant
#' number of components. Default is FALSE.
#' @return the VISION with latentSpace populated
computeLatentSpace <- function(object, projection_genes = NULL,
                               perm_wPCA = NULL) {

    message("Computing a latent space for expression data...\n")

    if (!is.null(projection_genes)) object@projection_genes <- projection_genes
    if (!is.null(perm_wPCA)) object@perm_wPCA <- perm_wPCA

    expr <- object@exprData
    weights <- object@weights
    projection_genes <- object@projection_genes
    perm_wPCA <- object@perm_wPCA

    if (!is.null(projection_genes)) {
        exprData <- expr[projection_genes, , drop = FALSE]
    } else {
        exprData <- expr
    }

    if (all(dim(weights) == c(1, 1))){
        weights <- matrix(1, nrow = nrow(exprData), ncol = ncol(exprData),
                          dimnames = list(
                                          rownames(exprData),
                                          colnames(exprData)
                                         )
                         )
    }

    exprData <- matLog2(exprData)

    if (perm_wPCA) {
        res <- applyPermutationWPCA(exprData, weights, components = 30)
        pca_res <- res[[1]]
    } else {
        res <- applyWeightedPCA(exprData, weights, maxComponents = 30)
        pca_res <- res[[1]]
    }

    object@latentSpace <- t(pca_res)
    return(object)
}

#' generate projections
#'
#' Generates 2-dimensional representations of the expression matrix
#' Populates the 'Projections' slot on the VISION object
#'
#'
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
generateProjections <- function(object) {
  message("Projecting data into 2 dimensions...")

  # Some projection methods operate on the full expression matrix
  # If using one of these, we need to compute 'projection_genes'
  projection_methods <- object@projection_methods
  if ("ICA" %in% projection_methods || "RBFPCA" %in% projection_methods) {
      object <- filterData(object)
  }

  projections <- generateProjectionsInner(object@exprData,
                                     object@latentSpace,
                                     projection_genes = object@projection_genes,
                                     projection_methods = object@projection_methods)

  # Add inputProjections
  for (proj in names(object@inputProjections)){
      projections[[proj]] <- object@inputProjections[[proj]]
  }

  object@Projections <- projections

  message("")

  return(object)
}

#' Compute local correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @param signatureBackground as returned by `calculateSignatureBackground`
#' @return the VISION object with values set for the analysis results
analyzeLocalCorrelations <- function(object, signatureBackground = NULL) {

  if (is.null(signatureBackground)){
      signatureBackground <- calculateSignatureBackground(object, num = 3000)
  }

  message("Evaluating local consistency of signatures in latent space...\n")

  sigConsistencyScores <- sigConsistencyScores(
                                object@latentSpace,
                                object@sigScores,
                                object@metaData,
                                signatureBackground)

  message("Clustering signatures...\n")
  sigClusters <- clusterSignatures(object@sigScores,
                                   object@metaData,
                                   sigConsistencyScores$fdr,
                                   sigConsistencyScores$sigProjMatrix,
                                   clusterMeta = object@pool)

  sigConsistencyScoresData <- ProjectionData(
                             Consistency = sigConsistencyScores$sigProjMatrix,
                             pValue = sigConsistencyScores$pVals,
                             FDR = sigConsistencyScores$fdr,
                             sigClusters = sigClusters)

  object@SigConsistencyScores <- sigConsistencyScoresData

  return(object)
}


#' Compute trajectory correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @param signatureBackground as returned by `calculateSignatureBackground`
#' @return the VISION object with values set for the analysis results
analyzeTrajectoryCorrelations <- function(object, signatureBackground = NULL) {

  if (is.null(signatureBackground)){
      signatureBackground <- calculateSignatureBackground(object, num = 3000)
  }

  message("Evaluating local consistency of signatures within trajectory model...\n")
  sigVTreeProj <- sigConsistencyScores(object@latentTrajectory,
                                       object@sigScores,
                                       object@metaData,
                                       signatureBackground)

  message("Clustering signatures...\n")
  sigTreeClusters <- clusterSignatures(object@sigScores,
                                       object@metaData,
                                       sigVTreeProj$pVals,
                                       sigVTreeProj$sigProjMatrix,
                                       clusterMeta = object@pool)

  TrajectoryConsistencyScores <- ProjectionData(
      Consistency = sigVTreeProj$sigProjMatrix,
      pValue = sigVTreeProj$pVals,
      FDR = sigVTreeProj$fdr,
      sigClusters = sigTreeClusters)

  object@TrajectoryConsistencyScores <- TrajectoryConsistencyScores

  return(object)
}


#' Compute Ranksums Test, for all factor meta data.  One level vs all others
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colRanks
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
clusterSigScores <- function(object) {

    message("Computing differential signature tests...\n")

    sigScores <- object@sigScores
    metaData <- object@metaData

    metaData <- metaData[rownames(sigScores), , drop = FALSE]

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

    if (ncol(sigScores) > 0){
        sigScoreRanks <- colRanks(sigScores,
                                  preserveShape = TRUE,
                                  ties.method = "average")
        dimnames(sigScoreRanks) <- dimnames(sigScores)
    } else {
        sigScoreRanks <- sigScores
    }

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

            r1 <- matrix_wilcox(sigScoreRanks, cluster_ii,
                                check_na = FALSE, check_ties = FALSE)

            r2 <- matrix_wilcox(numericMetaRanks, cluster_ii,
                                check_na = TRUE, check_ties = TRUE)

            r3 <- matrix_chisq(factorMeta, cluster_ii)

            pval <- c(r1$pval, r2$pval, r3$pval)
            stat <- c(r1$stat, r2$stat, r3$stat)
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

    object@ClusterSigScores <- out
    return(object)

}

#' Compute pearson correlation between signature scores and principle components
#'
#' Populations the PCAnnotatorData slot of the VISION object
#'
#' @importFrom pbmcapply pbmclapply
#' @param object the VISION object
#' @return pearsonCorr numeric matrix N_Signatures x N_PCs
calculatePearsonCorr <- function(object){

  message("Computing correlations between signatures and expression PCs...\n")
  sigMatrix <- object@sigScores
  metaData <- object@metaData
  latentSpace <- object@latentSpace

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
           if (is.na(pc_result$estimate)) {  # happens i std dev is 0 for a sig
               return(0)
           } else {
               return(pc_result$estimate)
           }
      })
      return(ls_col_cor)
  })

  pearsonCorr <- do.call(rbind, pearsonCorr)
  rownames(pearsonCorr) <- colnames(computedSigMatrix)

  pcaAnnotData <- PCAnnotatorData(pearsonCorr = pearsonCorr)
  object@PCAnnotatorData <- pcaAnnotData

  return(object)
}

convertToDense <- function(object) {

    object@exprData <- as.matrix(object@exprData)

    return(object)
}
