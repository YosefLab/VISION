#' create micro-clusters that reduce noise and complexity while maintaining
#' the overall signal in the data
#' @param object the FastProject object for which to cluster the cells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @param BPPARAM the parallelization backend to use
#' @return the FastProject with pooled cells
poolCells <- function(object,
                      cellsPerPartition=object@cellsPerPartition,
                      BPPARAM=SerialParam()) {
    object@cellsPerPartition <- cellsPerPartition
    microclusters <- applyMicroClustering(getExprData(object@exprData),
                                          object@housekeepingData,
                                          object@cellsPerPartition,
                                          BPPARAM)
    object@exprData <- ExpressionData(microclusters[[1]])
    object@pools <- microclusters[[2]]
    return(object)
}

#' filter data accourding to the provided filters
#' @param object the FastProject object
#' @param threshold threshold to apply for the threshold filter
#' @param filters list of filters to compute options are limited to
#' @return the FastProject object, populated with filtered data
filterData <- function(object,
                       threshold=object@threshold,
                       filters=object@filters) {
    object@filters <- filters
    object@threshold <- threshold
    message("filtering data...")
    if (object@threshold == 0) {
        num_samples <- ncol(getExprData(object@exprData))
        object@threshold <- round(0.2 * num_samples)
    }

    object@exprData <- applyFilters(object@exprData,
                                    object@threshold,
                                    object@filters)
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
    clustered <- (ncol(getExprData(object@exprData)) > 15000 || object@pool)
    if (!clustered && !object@nomodel) {
        message("Computing weights from False Negative Function...")
        falseneg_out <- createFalseNegativeMap(object@allData,
                                               object@housekeepingData)
        object@weights <- computeWeights(falseneg_out[[1]], falseneg_out[[2]],
                                         object@exprData)

    } else if (all(is.na(object@weights)) ||
               ncol(object@weights) != ncol(object@exprData)) {
        object@weights <- matrix(1L, nrow=nrow(object@allData),
                                 ncol=ncol(object@allData))
    }

    rownames(object@weights) <- rownames(object@allData)
    colnames(object@weights) <- colnames(object@allData)
    return(object)
}

#' Normalize the data
#' @param object the FastPrpject object
#' @param sig_norm_method (optional) Method to apply to normalize the expression
#' matrix before calculating signature scores
#' @return the FastProject object with normalized data matrices
normalizeData <- function(object,
                          sig_norm_method=object@sig_norm_method) {
    object@sig_norm_method <- sig_norm_method
    object@exprData <- updateExprData(object@exprData,
                                    getNormalizedCopy(object@exprData,
                                                      object@sig_norm_method))
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
#' @param precomputedData a list of existing cell-signatures
#' @param sig_score_method the scoring method to use
#' @param min_signature_genes the minimal number of genes that must be present
#' in both the data and the signature for the signature to be evaluated
#' @param BPPARAM the parallelization backend to use for these computations
#' @return the FastProject object, with signature score slots populated
calcSignatureScores <- function(object,
                                sigData=object@sigData,
                                precomputedData=object@precomputedData,
                                sig_score_method=object@sig_score_method,
                                min_signature_genes=object@min_signature_genes,
                                BPPARAM=NULL) {
    if(is.null(BPPARAM)) BPPARAM <- SerialParam()

    ## override object parameters
    if(!is.null(sigData)) object@sigData <- sigData
    if(!is.null(precomputedData)) object@precomputedData <- precomputedData
    object@sig_score_method <- sig_score_method
    object@min_signature_genes <- min_signature_genes

    sigScores <- bplapply(object@sigData, function(s) {
        singleSigEval(s, object@sig_score_method, object@exprData,
                      object@weights,
                      object@min_signature_genes)
    }, BPPARAM=BPPARAM)

    sigList <- object@sigData
    sigNames <- names(object@sigData)
    for (s in object@precomputedData) {
        if (length(s@sample_labels) != ncol(object@exprData)) {
            s@scores <- s@scores[colnames(object@exprData)]
            s@sample_labels <- colnames(object@exprData)
        }
        sigScores <- c(sigScores, s)
        sigList<- c(sigList, Signature(list(), s@name, "", "",
                                       isPrecomputed=TRUE,
                                       isFactor=s@isFactor, cluster=0))
        sigNames <- c(sigNames, s@name)
    }
    names(sigList) <- sigNames
    names(sigScores) <- sigNames

    # Remove any signatures that didn't compute correctly
    toRemove <- sapply(sigScores, is.null)
    sigScores <- sigScores[!toRemove]
    object@sigScores <- sigScores

    ## Convert Sig Scores to matrix
    names <- c()
    sigMatrix <- matrix(0L, nrow=length(sigScores),
                        ncol=length(sigScores[[1]]@scores))
    for (sig in 1:length(sigScores)) {
        names <- c(names, sigScores[[sig]]@name)
        scores <- t(as.matrix(sigScores[[sig]]@scores))
        sigMatrix[sig,] <- scores
    }

    rownames(sigMatrix) <- names
    colnames(sigMatrix) <- colnames(object@exprData)
    object@sigMatrix <- sigMatrix

    object@sigData <- sigList[rownames(sigMatrix)]

    # Construct random signatures for background distribution
    sigSizes <- lapply(object@sigData, function(s) length(s@sigDict))
    randomSigs <- generatePermutationNull(3000, object@exprData, sigSizes)

    ## Compute signature scores for random signatures generated
    randomSigScores <- bplapply(randomSigs, function(s) {
        singleSigEval(s, object@sig_score_method, object@exprData,
                      object@weights,
                      object@min_signature_genes)
    } ,BPPARAM=BPPARAM)
    names(randomSigScores) <- names(randomSigs)

    ## Remove random signatures that didn't compute correctly
    toRemove <- sapply(randomSigScores, is.null)
    object@randomSigScores <- randomSigScores[!toRemove]

    return(object)
}

#' Analyze projections
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
#' @param BPPARAM the parallelization backend to use
#' @return the FastProject object with values set for the analysis results
analyzeProjections <- function(object,
                               lean=object@lean,
                               perm_wPCA=object@perm_wPCA,
                               BPPARAM) {
  object@lean <- lean
  object@perm_wPCA <- perm_wPCA

  # Apply projections to filtered gene sets, create new projectionData object
  filterModuleList <- list()

  ## calculate projections - s
    for (filter in object@filters) {
        message("Filter: ", filter)

        message("Projecting data into 2 dimensions...")

        projectData <- generateProjections(object@exprData, object@weights,
                                           filter,
                                           inputProjections <- c(),
                                           lean=object@lean,
                                           perm_wPCA=object@perm_wPCA,
                                           BPPARAM = BPPARAM)

        message("Computing significance of signatures...")
        sigVProj <- sigsVsProjections(projectData$projections,
                                      object@sigScores,
                                      object@randomSigScores,
                                      BPPARAM=BPPARAM)

        message("Clustering Signatures...")
        sigClusters <- clusterSignatures(object@sigData, object@sigMatrix,
                                         sigVProj$pVals, k=10)

        projData <- ProjectionData(projections = projectData$projections,
                                   keys = sigVProj$projNames,
                                   sigProjMatrix = sigVProj$sigProjMatrix,
                                   pMatrix = sigVProj$pVals,
                                   sigClusters = sigClusters)

        message("Fitting principle tree...")
        treeProjs <- generateTreeProjections(object@exprData, filter,
                                   inputProjections = projectData$projections,
                                   permMats = projectData$permMats,
                                   BPPARAM = BPPARAM)

        message("Computing significance of signatures...")
        sigVTreeProj <- sigsVsProjections(treeProjs$projections,
                                          object@sigScores,
                                          object@randomSigScores,
                                          BPPARAM=BPPARAM)

        message("Clustering Signatures...")
        sigTreeClusters <- clusterSignatures(object@sigData, object@sigMatrix,
                                             sigVTreeProj$pVals, k=10)

        treeProjData <- TreeProjectionData(projections = treeProjs$projections,
                                   keys = sigVTreeProj$projNames,
                                   sigProjMatrix = sigVTreeProj$sigProjMatrix,
                                   pMatrix = sigVTreeProj$pVals,
                                   sigClusters = sigTreeClusters,
                                   treeScore = treeProjs$treeScore)

        ## remove the factors from the pearson correlation calculation
        precomp_names <- sapply(object@precomputedData, function(x) {
          ifelse(x@isFactor, x@name, NA)
        })
        precomp_names <- precomp_names[!is.na(precomp_names)]
        computedSigMatrix <- object@sigMatrix[setdiff(rownames(object@sigMatrix),
                                                      precomp_names),]


        pearsonCorr <- lapply(1:nrow(computedSigMatrix), function(i) {
          lapply(1:nrow(projectData$fullPCA), function(j) {
            return(calcPearsonCorrelation(computedSigMatrix[i,],
                                          projectData$fullPCA[j,]))
          })
        })

        pearsonCorr <- matrix(unlist(lapply(pearsonCorr, unlist)),
                              nrow=nrow(computedSigMatrix),
                              ncol=nrow(projectData$fullPCA), byrow=TRUE)

        rownames(pearsonCorr) <- rownames(computedSigMatrix)
        colnames(pearsonCorr) <- rownames(projectData$fullPCA)

        pcaAnnotData <- PCAnnotatorData(fullPCA = projectData$fullPCA,
                                        pearsonCorr = pearsonCorr,
                                        loadings = projectData$loadings)

        filterModuleData <- FilterModuleData(filter = filter,
                                             genes = projectData$geneNames,
                                             ProjectionData = projData,
                                             TreeProjectionData = treeProjData,
                                             PCAnnotatorData = pcaAnnotData)

        filterModuleList[[filter]] <- filterModuleData

    }

    object@filterModuleList <- filterModuleList
    return(object)
}
