#' Creates an initial clustering of the cells
#'
#' Results of this are stored as a new variabe in the object's metaData
#' and 'cluster_variable' is populated with its name
#'
#' @param object the FastProject object for which to cluster the cells
#' @return the FastProject object modifed as described above
clusterCells <- function(object) {

    exprData <- object@exprData
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

    n_workers <- getWorkerCount()

    kn <- ball_tree_knn(res,
                        min(round(sqrt(nrow(res))), 30),
                        n_workers)

    cl <- louvainCluster(kn, res)

    names(cl) <- paste('Cluster', seq(length(cl)))

    # cl is list of character vector
    cluster_variable <- "FastProject_Clusters"
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
#' @param object the FastProject object for which to cluster the cells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @return the FastProject with pooled cells
poolCells <- function(object,
                      cellsPerPartition=object@cellsPerPartition) {
    object@cellsPerPartition <- cellsPerPartition

    message(paste(
      "Performing micro-pooling on",
      ncol(object@exprData),
      "cells with a target pool size of",
      object@cellsPerPartition
    ))

    if (object@cluster_variable != "") {
        preserve_clusters <- object@metaData[[object@cluster_variable]]
        names(preserve_clusters) <- rownames(object@metaData)
    } else {
        preserve_clusters <- NULL
    }

    microclusters <- applyMicroClustering(object@exprData,
                                          cellsPerPartition = object@cellsPerPartition,
                                          filterInput = object@projection_genes,
                                          filterThreshold = object@threshold,
                                          preserve_clusters = preserve_clusters)

    object@exprData <- microclusters[[1]]
    object@pools <- microclusters[[2]]

    poolMeta <- createPooledMetaData(object@metaData, object@pools)

    object@metaData <- poolMeta

    return(object)
}

#' filter data accourding to the provided filters
#' @param object the FastProject object
#' @param threshold threshold to apply for the threshold filter
#' @param projection_genes either a list of genes or a method to select genes
#' @return the FastProject object, populated with filtered data
filterData <- function(object,
                       threshold=object@threshold,
                       projection_genes=object@projection_genes) {

    object@projection_genes <- projection_genes
    object@threshold <- threshold

    message("Determining Projection Genes...")

    if (object@threshold == 0) {
        num_samples <- ncol(object@exprData)
        object@threshold <- round(0.2 * num_samples)
    }

    object@projection_genes <- applyFilters(
                object@exprData,
                object@threshold,
                object@projection_genes)

    return(object)
}

#' Calculate weights based on the false negative rates, computed using the
#' provided housekeeping genes
#' @param object the FastProject object
#' @param nomodel (optional) if TRUE, no fnr curve calculated and all weights
#' equal to 1. Else FNR and weights calculated.
#' @return the FastProject object with populated weights slot
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

    } else if (all(is.na(object@weights)) ||
               ncol(object@weights) != ncol(object@exprData)) {
        object@weights <- matrix(1L, nrow=nrow(object@exprData),
                                 ncol=ncol(object@exprData))
    }

    rownames(object@weights) <- rownames(object@exprData)
    colnames(object@weights) <- colnames(object@exprData)
    return(object)
}

#' calculate signature scores
#'
#' For each signature-cell pair, compute a score that captures the level of
#' correspondence between the cell and the signature.
#' To estimate significance of these scores, a set of random gene signatures is
#' generated to create a null distribution
#'
#' @param object the FastProject object
#' @param sigData a list of Signature objects for which to compute the scores
#' @param metaData a list of existing cell-signatures
#' @param sig_norm_method (optional) Method to apply to normalize the expression
#' matrix before calculating signature scores
#' @param sig_score_method the scoring method to use
#' @param min_signature_genes the minimal number of genes that must be present
#' in both the data and the signature for the signature to be evaluated
#' @return the FastProject object, with signature score slots populated
calcSignatureScores <- function(object,
                                sigData=object@sigData,
                                metaData=object@metaData,
                                sig_norm_method=object@sig_norm_method,
                                sig_score_method=object@sig_score_method,
                                min_signature_genes=object@min_signature_genes) {

    message("Evaluating Signature Scores on Cells...")

    ## override object parameters
    if(!is.null(sigData)) object@sigData <- sigData
    if(!is.null(metaData)) object@metaData <- metaData
    object@sig_norm_method <- sig_norm_method
    object@sig_score_method <- sig_score_method
    object@min_signature_genes <- min_signature_genes

    normExpr <- getNormalizedCopy(object@exprData, object@sig_norm_method)

    sigScores <- batchSigEval(object@sigData, object@sig_score_method,
                              normExpr, object@weights, object@min_signature_genes)

    sigList <- object@sigData

    # Remove any signatures that didn't compute correctly
    toRemove <- vapply(sigScores, is.null, TRUE)
    sigScores <- sigScores[!toRemove]
    object@sigScores <- sigScores

    ## Convert Sig Scores to matrix
    # Note: This needs to be a data frame because some scores are factors
    object@sigMatrix <- data.frame(lapply(sigScores, function(x) x@scores),
                                  check.names = FALSE)

    object@sigData <- sigList[colnames(object@sigMatrix)]

    return(object)
}

#' analyze projections
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the FastProject object
#' @param lean if TRUE run a lean simulation. Else more robust pipeline
#' initiated. Default is FALSE
#' @param perm_wPCA If TRUE, apply permutation WPCA to calculate significant
#' number of PCs. Else not. Default FALSE.
#' @return the FastProject object with values set for the analysis results
analyzeProjections <- function(object,
                               lean=object@lean,
                               perm_wPCA=object@perm_wPCA) {

  object@lean <- lean
  object@perm_wPCA <- perm_wPCA

  message("Computing background distribution for signature scores...")
  signatureBackground <- calculateSignatureBackground(object, num=3000)

  message("Projecting data into 2 dimensions...")

  projectData <- generateProjections(object@exprData, object@weights,
                                     projection_genes=object@projection_genes,
                                     inputProjections=object@inputProjections,
                                     lean=object@lean,
                                     perm_wPCA=object@perm_wPCA)

  message("Evaluating signatures against projections...")
  sigVProj <- sigsVsProjections(projectData$projections,
                                object@sigScores,
                                object@metaData,
                                signatureBackground)

  message("Clustering Signatures...")
  sigClusters <- clusterSignatures(object@sigMatrix,
                                   object@metaData,
                                   sigVProj$pVals,
                                   clusterMeta = object@pool)

  projData <- ProjectionData(projections = projectData$projections,
                             keys = colnames(sigVProj$sigProjMatrix),
                             sigProjMatrix = sigVProj$sigProjMatrix,
                             pMatrix = sigVProj$pVals,
                             sigClusters = sigClusters,
                             emp_pMatrix = sigVProj$emp_pVals)

  if (tolower(object@trajectory_method) != "none") {
      message("Fitting principle tree...")
      treeProjs <- generateTreeProjections(projectData$fullPCA,
                                 inputProjections = projectData$projections,
                                 permMats = projectData$permMats)

      message("Computing significance of signatures...")
      sigVTreeProj <- sigsVsProjections(treeProjs$projections,
                                        object@sigScores,
                                        object@metaData,
                                        signatureBackground)

      message("Clustering Signatures...")
      sigTreeClusters <- clusterSignatures(object@sigMatrix,
                                           object@metaData,
                                           sigVTreeProj$pVals,
                                           clusterMeta = object@pool)

      treeProjData <- TreeProjectionData(projections = treeProjs$projections,
                                 keys = colnames(sigVTreeProj$sigProjMatrix),
                                 sigProjMatrix = sigVTreeProj$sigProjMatrix,
                                 pMatrix = sigVTreeProj$pVals,
                                 sigClusters = sigTreeClusters,
                                 treeScore = treeProjs$treeScore)
  } else {
      treeProjData <- NULL
  }

  message("Computing Correlations between Signatures and Expression PCs...")
  pearsonCorr <- calculatePearsonCorr(object@sigMatrix,
                     object@metaData, projectData$fullPCA)


  pcaAnnotData <- PCAnnotatorData(fullPCA = projectData$fullPCA,
                                  pearsonCorr = pearsonCorr,
                                  loadings = projectData$loadings)

  filterModuleData <- FilterModuleData(ProjectionData = projData,
                                       TreeProjectionData = treeProjData,
                                       PCAnnotatorData = pcaAnnotData)

  object@filterModuleData <- filterModuleData
  return(object)
}

#' Compute pearson correlation between signature scores and principle components
#'
#' @importFrom Hmisc rcorr
#' @param sigMatrix Signature scores dataframe
#' @param metaData data.frame of meta-data for cells
#' @param fullPCA numeric matrix N_Cells x N_PCs
#' @return pearsonCorr numeric matrix N_Signatures x N_PCs
calculatePearsonCorr <- function(sigMatrix, metaData, fullPCA){

  ## combined gene signs and numeric meta variables

  numericMetaVars <- vapply(colnames(metaData),
                            function(x) is.numeric(metaData[[x]]),
                            FUN.VALUE = TRUE)

  numericMeta <- metaData[, numericMetaVars]
  numericMeta <- numericMeta[rownames(sigMatrix), ]


  computedSigMatrix <- cbind(sigMatrix, numericMeta)

  pearsonCorr <- lapply(1:ncol(computedSigMatrix), function(i) {
    lapply(1:nrow(fullPCA), function(j) {
               ss <- computedSigMatrix[, i];
               pc <- fullPCA[j, ];
               pc_result <- rcorr(ss, pc, type="pearson")
               return(pc_result[[1]][1,2])
    })
  })

  pearsonCorr <- matrix(unlist(lapply(pearsonCorr, unlist)),
                        nrow=ncol(computedSigMatrix),
                        ncol=nrow(fullPCA), byrow=TRUE)

  rownames(pearsonCorr) <- colnames(computedSigMatrix)
  colnames(pearsonCorr) <- rownames(fullPCA)

  return(pearsonCorr)
}

convertToDense <- function(object) {

    object@exprData <- as.matrix(object@exprData)

    return(object)
}
