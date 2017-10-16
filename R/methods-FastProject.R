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
#' return(createUserGeneSignature(name = paste0("sig",i),
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
#' return(createUserGeneSignature(name = paste0("sig",i),
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

    ## compute signature scores - s
    # Score user defined signatures with defined method (default = weighted)
    sigScores <- c()
    if (object@sig_score_method == "naive") {
    message("Applying naive signature scoring method...")
    } else if (object@sig_score_method == "weighted_avg") {
    message("Applying weighted signature scoring method...")
    }

    sigScores <- bplapply(object@sigData, function(s) {
      singleSigEval(s, object@sig_score_method, object@exprData, object@weights,
                    object@min_signature_genes)
    }, BPPARAM=BPPARAM)

    sigSizes <- lapply(object@sigData, function(s) length(s@sigDict))

    sigList <- object@sigData
    sigNames <- names(object@sigData)
    for (s in object@precomputedData) {
        if (length(s@sample_labels) != ncol(object@exprData)) {
        s@scores <- s@scores[colnames(object@exprData)]
        s@sample_labels <- colnames(object@exprData)
    }
    sigScores <- c(sigScores, s)
    sigList<- c(sigList, Signature(list(), s@name, "", "", isPrecomputed=TRUE,
                                   isFactor=s@isFactor, cluster=0))
    sigNames <- c(sigNames, s@name)
    }
    names(sigList) <- sigNames
    names(sigScores) <- sigNames

    # Remove any signatures that didn't compute correctly
    toRemove <- sapply(sigScores, is.null)
    sigScores <- sigScores[!toRemove]

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


    # Construct random signatures for background distribution
    randomSigs <- generatePermutationNull(3000, object@exprData, sigSizes)


    ## Compute signature scores for random signatures generated
    randomSigScores <- c()
    if (object@sig_score_method == "naive") {
        message("Applying naive signature scoring method...")
    } else if (object@sig_score_method == "weighted_avg") {
        message("Applying weighted signature scoring method...")
    }

    randomSigScores <- bplapply(randomSigs, function(s) {
      singleSigEval(s, object@sig_score_method, object@exprData, object@weights,
                    object@min_signature_genes)
    } ,BPPARAM=BPPARAM)
    names(randomSigScores) <- names(randomSigs)

    ## Remove random signatures that didn't compute correctly
    toRemove <- sapply(randomSigScores, is.null)
    randomSigScores <- randomSigScores[!toRemove]

    ## calc sig scores - e
    # Apply projections to filtered gene sets, create new projectionData object
    filterModuleList <- list()

    ## calculate projections - s
    for (filter in object@filters) {
    message("Filter level: ", filter)

    message("Projecting data into 2 dimensions...")

    sigList <- sigList[rownames(sigMatrix)]

    projectData <- generateProjections(object@exprData, object@weights, filter,
                               inputProjections <- c(),
                               lean=object@lean, perm_wPCA = object@perm_wPCA,
                               BPPARAM = BPPARAM)
    projs <- projectData$projections ##TODO: not used
    g <- projectData$geneNames ## TODO: not used
    pca_res <- projectData$fullPCA
    loadings <- projectData$loadings ##TODO: not used
    permMats <- projectData$permMats ##TODO: not used

    message("Computing significance of signatures...")
    sigVProj <- sigsVsProjections(projectData$projections, sigScores,
                                  randomSigScores, BPPARAM=BPPARAM)

    message("Clustering Signatures...")
    sigClusters <- clusterSignatures(sigList, sigMatrix, sigVProj$pVals, k=10)

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
    sigVTreeProj <- sigsVsProjections(treeProjs$projections, sigScores,
                                        randomSigScores, BPPARAM=BPPARAM)

    message("Clustering Signatures...")
    sigTreeClusters <- clusterSignatures(sigList, sigMatrix,
                                         sigVTreeProj$pVals, k=10)

    treeProjData <- TreeProjectionData(projections = treeProjs$projections,
                                    keys = sigVTreeProj$projNames,
                                    sigProjMatrix = sigVTreeProj$sigProjMatrix,
                                    pMatrix = sigVTreeProj$pVals,
                                    sigClusters = sigTreeClusters,
                                    treeScore = treeProjs$treeScore)

    precomp_names <- sapply(object@precomputedData, function(x) ifelse(x@isFactor, x@name,))
    computedSigMatrix <- sigMatrix[setdiff(rownames(sigMatrix), precomp_names),]


    pearsonCorr <- lapply(1:nrow(computedSigMatrix), function(i) {
        lapply(1:nrow(pca_res), function(j) {
            return(calcPearsonCorrelation(computedSigMatrix[i,], pca_res[j,]))
            })
        })

        pearsonCorr <- lapply(pearsonCorr, unlist)
        pearsonCorr <- matrix(unlist(pearsonCorr), nrow=nrow(computedSigMatrix),
                            ncol=nrow(pca_res), byrow=TRUE)

        rownames(pearsonCorr) <- rownames(computedSigMatrix)
        colnames(pearsonCorr) <- rownames(pca_res)

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

    slots <- methods::slotNames("FastProject")
    fpParams <- list()
    for (s in slots) {
        fpParams[[s]] <- methods::slot(object, s)
    }

    ## TODO!!!
    # fpOut <- FastProjectOutput(eData, filterModuleList, sigMatrix, sigList,
    #                            fpParams, pools)
    #
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


normalizeData <- function(object) {
    normalizedData <- getNormalizedCopy(object@exprData, object@sig_norm_method)
    object@exprData <- updateExprData(object@exprData, normalizedData)
    return(object)
}

analyzeProjections <- function() {}
