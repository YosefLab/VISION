#' Initializes a new FastProject object.
#'
#' @import BiocParallel
#' @importFrom Biobase ExpressionSet exprs
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @import logging
#'
#' @param data expression data - can be one of these: \itemize{
#' \item numeric matrix
#' \item path of a file,
#' \item ExpressionSet object
#' \item SummerizedExperiment object (or extending classes)
#' }
#' @param signatures list of file paths to signature files (.gmt or .txt) or a
#' list of Signature objects
#' @param housekeeping file path to housekeeping data file, or vector of gene
#' names
#' @param norm_methods normalization methods to be extracted from the scone
#' object
#' @param precomputed data file with precomputed signature scores (.txt), or a
#' data.frame with meta information. Note that rows should match samples in the
#' data, and columns should either be factors or numerics.
#' @param nofilter if TRUE, no filter applied; else filters applied.
#' Default is FALSE
#' @param nomodel if TRUE, no fnr curve calculated and all weights equal to 1.
#' Else FNR and weights calculated.
#'              Default is TRUE.
#' @param filters list of filters to compute
#' @param lean if TRUE run a lean simulation. Else more robust pipeline
#' initiated. Default is FALSE
#' @param qc if TRUE calculate QC; else not. Default is FALSE
#' @param min_signature_genes Minimum number of genes required to compute a
#' signature
#' @param projections File containing precomputed projections for analysis
#' @param weights Precomputed weights for each coordinate. Normally computed
#' from the FNR curve.
#' @param threshold Threshold to apply for the threshold filter
#' @param perm_wPCA If TRUE, apply permutation WPCA to calculate significant
#' number of PCs. Else not. Default FALSE.
#' @param sig_norm_method Method to apply to normalize the expression matrix
#' before calculating signature scores
#' @param sig_score_method Method to apply when calculating signature scores
#' @param pool a boolean value indicating whether or not to create supercells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @return A FastProject object
#' @rdname FastProject-class
#' @export
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # create housekeeping genes
#' hkg <- paste0("gene",sample(1:500, 50))
#'
#' #create 20 signatures of 25 genes each
#' sigs <- lapply(1:20, function(i) {
#' sigData <- sign(rnorm(25))
#' names(sigData) <- paste0("gene",sample(1:100,25))
#' return(createGeneSignature(name = paste0("sig",i),
#'                                  sigData = sigData))
#' })
#'
#' fp <- FastProject(data = expMat,
#'                      signatures = sigs,
#'                      housekeeping = hkg)
setMethod("FastProject", signature(data = "matrix"),
            function(data, signatures, housekeeping=NULL, norm_methods = NULL,
                    precomputed=NULL, nofilter=FALSE, nomodel=FALSE,
                    filters=c("fano"), lean=FALSE, qc=FALSE,
                    min_signature_genes=5, projections="", weights=NULL,
                    threshold=0, perm_wPCA=FALSE, sig_norm_method="znorm_rows",
                    sig_score_method="weighted_avg", pool=FALSE,
                    cellsPerPartition=100) {

            .Object <- new("FastProject")

            rownames(data) <- sapply(rownames(data), toupper)
            .Object@allData = data
            .Object@exprData <- ExpressionData(data)

            if (is.null(housekeeping)) {
                .Object@housekeepingData <- c()
                .Object@nomodel = TRUE
            } else if (length(housekeeping) == 1) {
                .Object@housekeepingData <- readHKGToMatrix(housekeeping)
            } else {
                .Object@housekeepingData <- housekeeping
            }

            if (is.list(signatures)) {
                .Object@sigData <- signatures
            } else if (is.character(signatures)) {
                .Object@sigData <- readSignaturesInput(signatures)
            } else {
                stop("signatures must be paths to signature files or list of
                    Signature objects")
            }

            if (!is.null(precomputed)) {
                if(is.data.frame(precomputed)) {
                    .Object@precomputedData <-
                        SigScoresFromDataframe(precomputed)
                } else {
                    .Object@precomputedData <- readPrecomputed(
                    precomputed, colnames(.Object@exprData))
                }
            }

            if (is.null(weights)) {
                .Object@weights <- matrix(NA, nrow=10, ncol=0)
            } else {
                .Object@weights <- weights
            }

            .Object@nofilter <- nofilter
            if (!.Object@nomodel) {
                .Object@nomodel <- nomodel
            }
            .Object@filters <- filters
            .Object@projections <- projections
            .Object@threshold <- threshold
            .Object@sig_norm_method <- sig_norm_method
            .Object@sig_score_method <- sig_score_method
            .Object@lean = lean
            .Object@perm_wPCA = perm_wPCA
            .Object@pool = pool
            .Object@cellsPerPartition = cellsPerPartition

            return(.Object)
            }
)

#' @param ... additional arguments
#' @rdname FastProject-class
#' @export
setMethod("FastProject", signature(data = "character"),
            function(data, ...) {
            return(FastProject(readExprToMatrix(data), ...))
            }
)

#' @rdname FastProject-class
#' @export
setMethod("FastProject", signature(data = "ExpressionSet"),
            function(data, ...) {
            return(FastProject(Biobase::exprs(data), ...))
            }
)

#' @rdname FastProject-class
#' @export
setMethod("FastProject", signature(data = "SummarizedExperiment"),
          function(data, ...) {
            return(FastProject(SummarizedExperiment::assay(data), ...))
          }
)


#' Main entry point for running FastProject Analysis
#'
#' The main analysis function. Runs the entire FastProject analysis pipeline
#' and returns a FastProjectOutput object with the result,
#'
#' @export
#' @aliases Analyze
#' @param object FastProject object
#' @param BPPARAM a parallelization backend to use for the analysis
#' @return FastProjectOutput object
#'
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # create housekeeping genes
#' hkg <- paste0("gene",sample(1:500, 50))
#'
#' #create 20 signatures of 25 genes each
#' sigs <- lapply(1:20, function(i) {
#' sigData <- sign(rnorm(25))
#' names(sigData) <- paste0("gene",sample(1:100,25))
#' return(createGeneSignature(name = paste0("sig",i),
#'                                  sigData = sigData))
#' })
#'
#' fp <- FastProject(data = expMat,
#'                      housekeeping = hkg,
#'                      signatures = sigs)
#'
#' ## Analyze requires actual non-random data to run properly
#' \dontrun{
#' bp <- BiocParallel::SerialParam()
#' fp.out <- Analyze(fp, BPPARAM=bp)
#' }
setMethod("Analyze", signature(object="FastProject"),
            function(object, BPPARAM = NULL) {
    message("Beginning Analysis")
    if(is.null(BPPARAM)) {
        BPPARAM <- SerialParam()
    }

    object <- poolCells(object)

    object <- filterData(object)

    object <- calcWeights(object)

    object <- normalizeData(object)

    object <- calcSignatureScores(object, BPPARAM)

    object <- analyzeProjections(object, BPPARAM)

    slots <- methods::slotNames("FastProject")
    fpParams <- list()
    for (s in slots) {
        fpParams[[s]] <- methods::slot(object, s)
    }

    ## TODO!!!
    fpOut <- FastProjectOutput(filterModuleList,
                               fpParams)

    message("Analysis Complete!")

    #return(list(fpOut, timingList))
    return(fpOut)
})

#' create micro-clusters that reduce noise and complexity while maintaining
#' the overall signal in the data
#' @param object the FastProject object to cluster the cells of
#' @return the FastProject with pooled cells
poolCells <- function(object) {
    if (ncol(getExprData(object@exprData)) > 15000 || object@pool) {
        microclusters <- applyMicroClustering(getExprData(object@exprData),
                                              object@housekeepingData,
                                              object@cellsPerPartition,
                                              BPPARAM)
        object@exprData <- ExpressionData(microclusters[[1]])
        object@pools <- microclusters[[2]]
    }
    return(object)
}

#' filter data accourding to the provided filters
#' @param object the FastProject object
#' @return the FastProject object, populated with filtered data
filterData <- function(object) {
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
#' @return the FastProject object with populated weights slot
calcWeights <- function(object) {
    ##TODO:patch
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
#' @return the FastProject object with normalized data matrices
normalizeData <- function(object) {
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
#' @param BPPARAM the parallelization backend to use for these computations
#' @return the FastProject object, with signature score slots populated
calcSignatureScores <- function(object, BPPARAM=NULL) {
    if(is.null(BPPARAM)) BPPARAM <- SerialParam()

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
    object@sigScores <- sigScores[!toRemove]

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
      singleSigEval(s, object@sig_score_method, object@exprData, object@weights,
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
#' @param BPPARAM the parallelization backend to use
#' @return the FastProject object with values set for the analysis results
analyzeProjections <- function(object, BPPARAM) {
  # Apply projections to filtered gene sets, create new projectionData object
  filterModuleList <- list()

  ## calculate projections - s
  for (filter in object@filters) {
    message("Filter: ", filter)

    message("Projecting data into 2 dimensions...")

    projectData <- generateProjections(object@exprData, object@weights, filter,
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
    sigVTreeProj <- sigsVsProjections(treeProjs$projections, object@sigScores,
                                      object@randomSigScores, BPPARAM=BPPARAM)

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
