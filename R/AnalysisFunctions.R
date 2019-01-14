#' Creates an initial clustering of the cells
#'
#' Results of this are stored as a new variabe in the object's metaData
#' and 'cluster_variable' is populated with its name
#'
#' @param object the VISION object for which to cluster the cells
#' @return the VISION object modifed as described above
clusterCells <- function(object) {

    message("Clustering cells...")

    if (sum(dim(object@latentSpace)) == 2) { # No latent Space

        exprData <- matLog2(object@exprData)

        filterInput <- object@projection_genes
        filterThreshold <- object@threshold

        gene_passes <- applyFilters(exprData, filterThreshold, filterInput)
        fexpr <- exprData[gene_passes, ]

        # Compute wcov using matrix operations to avoid
        # creating a large dense matrix

        N <- ncol(fexpr)
        wcov <- tcrossprod(fexpr) / N

        mu <- as.matrix(rowMeans(fexpr), ncol = 1)
        mumu <- tcrossprod(mu)
        wcov <- as.matrix(wcov - mumu)

        # SVD of wieghted correlation matrix
        ncomp <- min(ncol(fexpr), nrow(fexpr), 10)
        decomp <- rsvd::rsvd(wcov, k = ncomp)
        evec <- t(decomp$u)

        # Project down using computed eigenvectors
        res <- (evec %*% fexpr) - as.vector(evec %*% mu)
        res <- as.matrix(res)
        res <- t(res) # avoid transposing many times below
    } else {
        res <- object@latentSpace
    }

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
    object@cluster_variable <- cluster_variable

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
      message("Performing micro-poolong on pre-computed pools")
    }

    if (object@cluster_variable != "") {
        preserve_clusters <- object@metaData[[object@cluster_variable]]
        names(preserve_clusters) <- rownames(object@metaData)
    } else {
        preserve_clusters <- NULL
    }

    if (length(object@pools) == 0) {
        pools <- applyMicroClustering(object@exprData,
                                              cellsPerPartition = object@cellsPerPartition,
                                              filterInput = object@projection_genes,
                                              filterThreshold = object@threshold,
                                              preserve_clusters = preserve_clusters,
                                              latentSpace = object@latentSpace)

        object@pools <- pools
    }

    pooled_cells <- poolMatrixCols(object@exprData, object@pools)
    object@exprData <- pooled_cells

    pooled_unnorm <- poolMatrixCols(object@unnormalizedData, object@pools)
    object@unnormalizedData <- pooled_unnorm

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

    message("Determining Projection Genes...")

    exprData <- matLog2(object@exprData)

    object@projection_genes <- applyFilters(
                exprData,
                object@threshold,
                object@projection_genes)

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

    message("Evaluating Signature Scores on Cells...")

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

    message("Computing a latent space for expression data...")

    if (!is.null(projection_genes)) object@projection_genes <- projection_genes
    if (!is.null(perm_wPCA)) object@perm_wPCA <- perm_wPCA

    expr <- object@exprData
    weights <- object@weights
    projection_genes <- object@projection_genes
    perm_wPCA <- object@perm_wPCA

    if (!is.null(projection_genes)) {
        exprData <- expr[projection_genes, ]
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

  projections <- generateProjectionsInner(object@exprData,
                                     object@latentSpace,
                                     projection_genes = object@projection_genes,
                                     projection_methods = object@projection_methods)

  # Add inputProjections
  for (proj in names(object@inputProjections)){
      projections[[proj]] <- object@inputProjections[[proj]]
  }

  object@Projections <- projections

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
      message("Computing background distribution for signature scores...")
      signatureBackground <- calculateSignatureBackground(object, num = 3000)
  }

  message("Evaluating local consistency of signatures...")

  sigConsistencyScores <- sigConsistencyScores(
                                object@latentSpace,
                                object@sigScores,
                                object@metaData,
                                signatureBackground)

  message("Clustering Signatures...")
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
      message("Computing background distribution for signature scores...")
      signatureBackground <- calculateSignatureBackground(object, num = 3000)
  }

  message("Computing significance of signatures...")
  sigVTreeProj <- sigConsistencyScores(object@latentTrajectory,
                                       object@sigScores,
                                       object@metaData,
                                       signatureBackground)

  message("Clustering Signatures...")
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
#' @importFrom parallel mclapply
#' @importFrom matrixStats colRanks
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
clusterSigScores <- function(object) {

    message("Computing differential signature tests...")

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

    out <- lapply(clusterMeta, function(variable){
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
    })

    object@ClusterSigScores <- out
    return(object)

}

#' Compute pearson correlation between signature scores and principle components
#'
#' Populations the PCAnnotatorData slot of the VISION object
#'
#' @param object the VISION object
#' @return pearsonCorr numeric matrix N_Signatures x N_PCs
calculatePearsonCorr <- function(object){

  message("Computing Correlations between Signatures and Expression PCs...")
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

  pearsonCorr <- matrix(
      0,
      nrow = ncol(computedSigMatrix),
      ncol = ncol(latentSpace),
      dimnames = list(
          colnames(computedSigMatrix),
          colnames(latentSpace)
      )
  )

  for (i in seq_len(ncol(computedSigMatrix))) {
      for (j in seq_len(ncol(latentSpace))) {
           ss <- computedSigMatrix[, i];
           pc <- latentSpace[, j];
           suppressWarnings({
               pc_result <- cor.test(ss, pc)
           })
           if (is.na(pc_result$estimate)) {  # happens is std dev is 0 for a sig
               pearsonCorr[i, j] <- 0
           } else {
               pearsonCorr[i, j] <- pc_result$estimate
           }
      }
  }

  pcaAnnotData <- PCAnnotatorData(pearsonCorr = pearsonCorr)
  object@PCAnnotatorData <- pcaAnnotData

  return(object)
}

convertToDense <- function(object) {

    object@exprData <- as.matrix(object@exprData)

    return(object)
}
