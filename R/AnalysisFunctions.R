#' create micro-clusters that reduce noise and complexity while maintaining
#' the overall signal in the data
#' @param object the FastProject object for which to cluster the cells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @return the FastProject with pooled cells
poolCells <- function(object,
                      cellsPerPartition=object@cellsPerPartition) {
    object@cellsPerPartition <- cellsPerPartition
    microclusters <- applyMicroClustering(getExprData(object@exprData),
                                          hkg=object@housekeepingData,
                                          cellsPerPartition=object@cellsPerPartition)
    object@exprData <- ExpressionData(microclusters[[1]])
    object@pools <- microclusters[[2]]
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
        num_samples <- ncol(getExprData(object@exprData))
        object@threshold <- round(0.2 * num_samples)
    }

    object@exprData <- applyFilters(object@exprData,
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
    clustered <- (ncol(getExprData(object@exprData)) > 15000 || object@pool)
    if (!clustered && !object@nomodel) {
        message("Computing weights from False Negative Rate curves...")
        falseneg_out <- createFalseNegativeMap(getExprData(object@exprData),
                                               object@housekeepingData)
        object@weights <- computeWeights(falseneg_out[[1]], falseneg_out[[2]],
                                         object@exprData)

    } else if (all(is.na(object@weights)) ||
               ncol(object@weights) != ncol(object@exprData)) {
        object@weights <- matrix(1L, nrow=nrow(getExprData(object@exprData)),
                                 ncol=ncol(getExprData(object@exprData)))
    }

    rownames(object@weights) <- rownames(getExprData(object@exprData))
    colnames(object@weights) <- colnames(getExprData(object@exprData))
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
#' @param sig_norm_method (optional) Method to apply to normalize the expression
#' matrix before calculating signature scores
#' @param sig_score_method the scoring method to use
#' @param min_signature_genes the minimal number of genes that must be present
#' in both the data and the signature for the signature to be evaluated
#' @return the FastProject object, with signature score slots populated
calcSignatureScores <- function(object,
                                sigData=object@sigData,
                                precomputedData=object@precomputedData,
                                sig_norm_method=object@sig_norm_method,
                                sig_score_method=object@sig_score_method,
                                min_signature_genes=object@min_signature_genes) {

    message("Evaluating Signature Scores on Cells...")

    ## override object parameters
    if(!is.null(sigData)) object@sigData <- sigData
    if(!is.null(precomputedData)) object@precomputedData <- precomputedData
    object@sig_norm_method <- sig_norm_method
    object@sig_score_method <- sig_score_method
    object@min_signature_genes <- min_signature_genes

    normExpr <- getNormalizedCopy(object@exprData, object@sig_norm_method)

    sigScores <- batchSigEval(object@sigData, object@sig_score_method,
                              normExpr, object@weights, object@min_signature_genes)

    sigList <- object@sigData
    for (s in object@precomputedData) {

        if(is.null(object@pools)){
            sigScores[[s@name]] <- s
        } else {
            if(s@isFactor){
                ## Need to compute majority level in each group
                N_MicroClusters <- ncol(normExpr)
                newlevels <- union("~", levels(s@scores))

                clustScores <- factor(integer(N_MicroClusters),
                                levels=newlevels)
                names(clustScores) <- colnames(normExpr)
                for(clust in names(clustScores)){

                    pool <- object@pools[[clust]]

                    if(length(pool) == 0){ #TODO: This shouldn't happen
                        clustScores[clust] <- "~"
                        next
                    }

                    vals <- s@scores[match(pool, names(s@scores))]
                    freq <- table(vals) / length(vals)
                    maxval <- freq[which.max(freq)]
                    if(maxval >= .5){
                        clust_val <- names(maxval)
                    } else {
                        clust_val <- "~"
                    }
                    clustScores[clust] <- clust_val
                }
            } else { # Then it must be numeric, just average

                ## Need to compute majority level in each group
                N_MicroClusters <- ncol(normExpr)

                clustScores <- numeric(N_MicroClusters)
                names(clustScores) <- colnames(normExpr)
                for(clust in names(clustScores)){

                    pool <- object@pools[[clust]]

                    if(length(pool) == 0){ #TODO: This shouldn't happen
                        clustScores[clust] = 0
                        next
                    }

                    vals <- s@scores[match(pool, names(s@scores))]
                    clust_val = mean(vals)
                    clustScores[clust] = clust_val
                }
            }

            sigScores[[s@name]] <- SignatureScores(
                                       clustScores, s@name, s@isFactor,
                                       s@isPrecomputed, 0)
        }


        sigList[[s@name]] <- Signature(list(), s@name, "", "",
                                       isPrecomputed=TRUE,
                                       isFactor=s@isFactor, cluster=0)
    }

    # Remove any signatures that didn't compute correctly
    toRemove <- vapply(sigScores, is.null, TRUE)
    sigScores <- sigScores[!toRemove]
    object@sigScores <- sigScores

    ## Convert Sig Scores to matrix
    # Note: This needs to be a data frame because some scores are factors
    object@sigMatrix = data.frame(lapply(sigScores, function(x) x@scores),
                                  check.names=FALSE)

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
  randomSigScores <- calculateSignatureBackground(object, num=3000)

  # Apply projections to filtered gene sets, create new projectionData object
  filterModuleList <- list()
  filter <- "fano"

  message("Projecting data into 2 dimensions...")

  projectData <- generateProjections(object@exprData, object@weights,
                                     filter,
                                     inputProjections=object@inputProjections,
                                     lean=object@lean,
                                     perm_wPCA=object@perm_wPCA)

  message("Evaluating signatures against projections...")
  sigVProj <- sigsVsProjections(projectData$projections,
                                object@sigScores,
                                randomSigScores)

  message("Clustering Signatures...")
  sigClusters <- clusterSignatures(object@sigData, object@sigMatrix,
                                   sigVProj$pVals)

  projData <- ProjectionData(projections = projectData$projections,
                             keys = colnames(sigVProj$sigProjMatrix),
                             sigProjMatrix = sigVProj$sigProjMatrix,
                             pMatrix = sigVProj$pVals,
                             sigClusters = sigClusters,
                             emp_pMatrix = sigVProj$emp_pVals)

  if (tolower(object@trajectory_method) != "none") {
      message("Fitting principle tree...")
      treeProjs <- generateTreeProjections(projectData$fullPCA, filter,
                                 inputProjections = projectData$projections,
                                 permMats = projectData$permMats)

      message("Computing significance of signatures...")
      sigVTreeProj <- sigsVsProjections(treeProjs$projections,
                                        object@sigScores,
                                        randomSigScores)

      message("Clustering Signatures...")
      sigTreeClusters <- clusterSignatures(object@sigData, object@sigMatrix,
                                           sigVTreeProj$pVals)

      treeProjData <- TreeProjectionData(projections = treeProjs$projections,
                                 keys = colnames(sigVTreeProj$sigProjMatrix),
                                 sigProjMatrix = sigVTreeProj$sigProjMatrix,
                                 pMatrix = sigVTreeProj$pVals,
                                 sigClusters = sigTreeClusters,
                                 treeScore = treeProjs$treeScore)
  } else {
      treeProjData <- NULL
  }

  ## remove the factors from the pearson correlation calculation
  factor_precomp <- vapply(object@precomputedData,
                           function(x) x@isFactor, FUN.VALUE=TRUE)
  precomp_names <- vapply(object@precomputedData,
                           function(x) x@name, FUN.VALUE="")

  precomp_names <- precomp_names[factor_precomp]

  computedSigMatrix <- object@sigMatrix[,setdiff(colnames(object@sigMatrix),
                                                precomp_names),drop=FALSE]
  computedSigMatrix <- t(as.matrix(computedSigMatrix))


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


  object@filterModuleList <- filterModuleList
  return(object)
}

convertToDense <- function(object) {

    object@exprData@data <- as.matrix(object@exprData@data)
    object@exprData@fanoFilter <- as.matrix(object@exprData@fanoFilter)

    return(object)
}
