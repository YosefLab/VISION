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
                        precomputed, colnames(.Object@allData))
                } else {
                    .Object@precomputedData <- readPrecomputed(
                        precomputed, colnames(.Object@allData))
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
#' @return FastProject object
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

    message("Analysis Complete!")

    return(object)
})

#' Save the FastProjectOutput object as an .RDS file and view the results on a
#' localhost
#'
#' Save the results object as an RDS file for future use, and launch a local
#' server to explore the results with a browser.
#'
#' @param fpout FastProject object
#' @param ofile the path to save the object in. If NULL, the object is saved
#' in the working directory [default:NULL]
#' @return the path of the saved file
#' @aliases saveFPOutAndViewResults
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
#'                      housekeeping = hkg,
#'                      signatures = sigs)
#'
#' ## Analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- Analyze(fp)
#' saveFPOutAndViewResults(fp.out)
#' }
setMethod("saveFPOutAndViewResults", signature(fpout="FastProject"),
          function(fpout, ofile=NULL) {
            if(is.null(ofile)) {
              i <- 1
              ofile <- paste0("./fpout", i, ".rds")
              while (file.exists(ofile)) {
                i <- i+1
                ofile <- paste0("./fpout", i, ".rds")
              }
            }

            saveRDS(fpout, file=ofile)
            viewResults(fpout)
            return(ofile)
          })

#' View results of analysis without saving output object
#'
#' launch a local server to explore the results with a browser.
#'
#' @param object FastProject object or path to a file containing such an
#' object (saved using saveAndViewResults, or directly using saveRDS)
#' @aliases viewResults
#' @return None
#' @export
#' @rdname viewResults
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
#' fp.out <- Analyze(fp)
#' viewResults(fp.out)
#' }
setMethod("viewResults", signature(object="FastProject"),
          function(object) {

            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            path <- find.package("FastProjectR")
            curpath <- getwd()
            tryCatch(expr = {
              e <- new.env()
              e$arg1 <- object
              sys.source(file = file.path(path, "FastProjectR_Output",
                                          "server_script.R"), envir = e)
            }, finally = {
              setwd(curpath)
            })
          })

#' @rdname viewResults
#' @export
setMethod("viewResults", signature(object="character"),
          function(object) {
            fpo <- readRDS(object)
            if(!methods::is(fpo, "FastProject")){
              stop("loaded object not a valid FastProject object")
            }
            viewResults(fpo)
          })

