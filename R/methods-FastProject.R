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
            .Object@exprData <- data
            rownames(.Object@exprData) <- sapply(rownames(.Object@exprData),
                                                 toupper)

            if (is.null(housekeeping)) {
                .Object@housekeepingData <- character()
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
                    .Object@precomputedData <- SigScoresFromDataframe(
                        precomputed, colnames(.Object@exprData))
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
            .Object@allData = .Object@exprData
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

    ptm <- Sys.time()
    timingList <- (ptm - ptm)
    tRows <- c("Start")

    clustered <- FALSE
    pools <- list()
    if (ncol(object@exprData) > 15000 || object@pool) {
        microclusters <- applyMicroClustering(object@exprData,
                                            object@housekeepingData,
                                            object@cellsPerPartition,
                                            BPPARAM)
        pooled_cells <- microclusters[[1]]
        pools <- microclusters[[2]]

        object@exprData <- pooled_cells
        clustered <- TRUE
    }

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Partition & Sample")

    # Wrap expression data frame into a ExpressionData class
    eData <- ExpressionData(object@exprData)

    # If no filter threshold was specified, set it to 20% of samples
    if (object@threshold == 0) {
    num_samples <- ncol(getExprData(eData))
    object@threshold <- round(0.2 * num_samples)
    }

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Threshold")

    filtered <- applyFilters(eData, object@threshold, object@filters)
    eData <- filtered[[1]]
    filterList <- filtered[[2]]

    originalData <- getExprData(eData)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Filter")

    if (!clustered && !object@nomodel) {
    falseneg_out <- createFalseNegativeMap(originalData, object@housekeepingData)
    func <- falseneg_out[[1]]
    params <- falseneg_out[[2]]

    message("Computing weights from False Negative Function...")
    object@weights <- computeWeights(func, params, eData)

    } else if (all(is.na(object@weights)) || ncol(object@weights) != ncol(object@exprData)) {
    object@weights <- matrix(1L, nrow=nrow(originalData), ncol=ncol(originalData))
    }

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "FN Function")

    rownames(object@weights) <- rownames(originalData)
    colnames(object@weights) <- colnames(originalData)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Weights")

    normalizedData <- getNormalizedCopy(eData, object@sig_norm_method)
    eData <- updateExprData(eData, normalizedData)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Normalize")

    # Score user defined signatures with defined method (default = weighted)
    sigScores <- c()
    if (object@sig_score_method == "naive") {
    message("Applying naive signature scoring method...")
    } else if (object@sig_score_method == "weighted_avg") {
    message("Applying weighted signature scoring method...")
    }

    sigScores <- bplapply(object@sigData, function(s) {
      singleSigEval(s, object@sig_score_method, eData, object@weights,
                    object@min_signature_genes)
    }, BPPARAM=BPPARAM)

    sigSizes <- lapply(object@sigData, function(s) length(s@sigDict))

    sigList <- object@sigData
    sigNames <- names(object@sigData)
    ## TODO: get rid of this somehow
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
    names <- c(names, toupper(sigScores[[sig]]@name))
    scores <- t(as.matrix(sigScores[[sig]]@scores))
    sigMatrix[sig,] <- scores
    }

    rownames(sigMatrix) <- names
    colnames(sigMatrix) <- colnames(object@exprData)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Sig Scores")

    # Construct random signatures for background distribution
    randomSigs <- generatePermutationNull(3000, eData, sigSizes)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Random Sigs")

    ## Compute signature scores for random signatures generated
    randomSigScores <- c()
    if (object@sig_score_method == "naive") {
    message("Applying naive signature scoring method...")
    } else if (object@sig_score_method == "weighted_avg") {
    message("Applying weighted signature scoring method...")
    }

    randomSigScores <- bplapply(randomSigs, function(s) {
      singleSigEval(s, object@sig_score_method, eData, object@weights,
                    object@min_signature_genes)
    } ,BPPARAM=BPPARAM)
    names(randomSigScores) <- names(randomSigs)

    ## Remove random signatures that didn't compute correctly
    toRemove <- sapply(randomSigScores, is.null)
    randomSigScores <- randomSigScores[!toRemove]

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Rand Sig Scores")

    # Apply projections to filtered gene sets, create new projectionData object
    filterModuleList <- list()
    for (filter in filterList) {
    message("Filter level: ", filter)

    message("Projecting data into 2 dimensions...")

    sigList <- sigList[rownames(sigMatrix)]

    projectData <- generateProjections(eData, object@weights, filter,
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

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("Pr. ", filter))

    message("Clustering Signatures...")
    sigClusters <- clusterSignatures(sigList, sigMatrix, sigVProj$pVals, k=10)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("ClusterSignaturesProjections ", filter))

    projData <- ProjectionData(projections = projectData$projections,
                                keys = sigVProj$projNames,
                                sigProjMatrix = sigVProj$sigProjMatrix,
                                pMatrix = sigVProj$pVals,
                                sigClusters = sigClusters)

    message("Fitting principle tree...")
    treeProjs <- generateTreeProjections(eData, filter,
                                    inputProjections = projectData$projections,
                                    permMats = projectData$permMats,
                                    BPPARAM = BPPARAM)
    message("Computing significance of signatures...")
    sigVTreeProj <- sigsVsProjections(treeProjs$projections, sigScores,
                                        randomSigScores, BPPARAM=BPPARAM)

    message("Clustering Signatures...")
    sigTreeClusters <- clusterSignatures(sigList, sigMatrix,
                                         sigVTreeProj$pVals, k=10)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("ClusterSignaturesProjections ", filter))

    treeProjData <- TreeProjectionData(projections = treeProjs$projections,
                                    keys = sigVTreeProj$projNames,
                                    sigProjMatrix = sigVTreeProj$sigProjMatrix,
                                    pMatrix = sigVTreeProj$pVals,
                                    sigClusters = sigTreeClusters,
                                    treeScore = treeProjs$treeScore)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("SigVTreePr ", filter))

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

        timingList <- rbind(timingList, c(difftime(Sys.time(), ptm,
                                                   units="secs")))
        tRows <- c(tRows, paste("Pearson Correlation", filter))

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

    fpOut <- FastProjectOutput(eData, filterModuleList, sigMatrix, sigList,
                               fpParams, pools)

    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, "Final")

    timingList <- as.matrix(timingList)
    rownames(timingList) <- tRows

    message("Analysis Complete!")

    #return(list(fpOut, timingList))
    return(fpOut)
})
