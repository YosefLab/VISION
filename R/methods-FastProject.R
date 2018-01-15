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
#' @param housekeeping vector of gene names
#' @param norm_methods normalization methods to be extracted from the scone
#' object
#' @param precomputed data file with precomputed signature scores (.txt), or a
#' data.frame with meta information. Note that rows should match samples in the
#' data, and columns should either be factors or numerics.
#' @param nomodel if TRUE, no fnr curve calculated and all weights equal to 1.
#' Else FNR and weights calculated. [Default:FALSE]
#' @param projection_genes name of method ('threshold' or 'fano') or list of 
#' genes to use when computing projections.
#' @param lean if TRUE run a lean simulation. Else more robust pipeline
#' initiated. Default is FALSE
#' @param min_signature_genes Minimum number of genes required to compute a
#' signature
#' @param weights Precomputed weights for each coordinate. Normally computed
#' from the FNR curve.
#' @param threshold Threshold to apply for the threshold filter
#' @param perm_wPCA If TRUE, apply permutation WPCA to calculate significant
#' number of PCs. Else not. Default FALSE.
#' @param sig_norm_method Method to apply to normalize the expression matrix
#' before calculating signature scores
#' @param sig_score_method Method to apply when calculating signature scores
#' @param pool indicates whether or not to create supercells. Acceptable values
#' are TRUE, FALSE, or 'auto', the last of which is the default and enables
#' pooling if there are more than 15000 cells.
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @param name a name for the sample - shown on the output report
#' @return A FastProject object
#' @rdname FastProject-class
#' @export
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
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
setMethod("FastProject", signature(data = "matrixORSparse"),
            function(data, signatures, housekeeping=NULL, norm_methods = NULL,
                    precomputed=NULL, nomodel=FALSE,
                    projection_genes=c("fano"), lean=FALSE, min_signature_genes=5,
                    weights=NULL, threshold=0, perm_wPCA=FALSE,
                    sig_norm_method="znorm_rows",
                    sig_score_method="weighted_avg", pool="auto",
                    cellsPerPartition=100, name=NULL) {

            .Object <- new("FastProject")

            rownames(data) <- toupper(rownames(data))
            .Object@allData = data
            .Object@exprData <- ExpressionData(data)

            if (is.null(housekeeping)) {
                .Object@housekeepingData <- character()
                .Object@nomodel = TRUE
            } else {
                .Object@housekeepingData <- housekeeping
            }

            if (is.list(signatures)) {
                .Object@sigData <- signatures
                names(.Object@sigData) <- vapply(.Object@sigData,
                                            function(x){x@name}, "")
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

            if (!.Object@nomodel) {
                .Object@nomodel <- nomodel
            }
            .Object@projection_genes <- projection_genes
            .Object@threshold <- threshold
            .Object@sig_norm_method <- sig_norm_method
            .Object@sig_score_method <- sig_score_method
            .Object@lean = lean
            .Object@perm_wPCA = perm_wPCA

            LOTS_OF_CELLS <- ncol(getExprData(.Object@exprData)) > 15000

            if (is.character(pool))
            {
                pool <- tolower(pool)
                if (pool == 'auto')
                {
                    if(LOTS_OF_CELLS) {
                        pool <- TRUE
                        message(paste(
                          "Over 15000 input cells detected.  Enabling micropooling with max",
                          cellsPerPartition, "cells per pool."
                          )
                        )
                    } else {
                        pool <- FALSE
                    }
                } else {
                    stop("Bad value for 'pool' argument")
                }
            }

            if (LOTS_OF_CELLS && !pool) {
                message(paste(
                      "Warning: Input data consists of",
                      ncol(getExprData(.Object@exprData)),
                      "cells and pool=FALSE.  It is recommend to set pool=TRUE when running on large numbers of cells"
                      )
                )
            }

            .Object@pool = pool
            .Object@cellsPerPartition = cellsPerPartition

            if (!is.null(name)) {
                .Object@name <- name
            }

            return(.Object)
            }
)

#' @param ... additional arguments
#' @rdname FastProject-class
#' @export
setMethod("FastProject", signature(data = "character"),
            function(data, ...) {
            return(FastProject(readExprAsMatrix(data), ...))
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
#' and returns a FastProject object populated with the result,
#'
#' @export
#' @aliases analyze
#' @param object FastProject object
#' @param BPPARAM a parallelization backend to use for the analysis
#' @return FastProject object
#'
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
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
#' ## analyze requires actual non-random data to run properly
#' \dontrun{
#' bp <- BiocParallel::SerialParam()
#' fp.out <- analyze(fp, BPPARAM=bp)
#' }
setMethod("analyze", signature(object="FastProject"),
            function(object, BPPARAM = NULL) {
    message("Beginning Analysis")
    if(is.null(BPPARAM)) {
        BPPARAM <- SerialParam()
    }

    if (object@pool) {
        object <- poolCells(object, BPPARAM = BPPARAM)
    }

    object <- filterData(object)
    object <- convertToDense(object)

    object <- calcWeights(object)

    object <- calcSignatureScores(object, BPPARAM = BPPARAM)

    object <- analyzeProjections(object, BPPARAM = BPPARAM)

    message("Analysis Complete!")

    return(object)
})

#' Add a set of projection coordinates
#'
#' @export
#' @param object FastProject object
#' @param name Name of the projection
#' @param coordinates numeric matrix or data.frame. Coordinates of each
#' sample in the projection (NUM_SAMPLES x NUM_COMPONENTS)
#' @return FastProject object
setMethod("addProjection", signature(object="FastProject"),
            function(object, name, coordinates) {

    if(is(coordinates, "data.frame")){
        coordinates = as.matrix(coordinates)
    }

    # Verify that projection coordinates are correct
    samples = object@exprData@data
    sample_names = colnames(samples)

    if(length(intersect(sample_names, rownames(coordinates))) != dim(coordinates)[1]){
        stop("Supplied coordinates must have rowlabels that match sample/cell names")
    }

    if(dim(coordinates)[2] != 2){
        stop("Projection must have exactly 2 components")
    }


    # Add it to the object
    proj = Projection(name, coordinates)

    object@inputProjections = c(object@inputProjections, proj)

    return(object)
})

#' Save the FastProject object as an .RDS file and view the results on a
#' localhost
#'
#' Save the results object as an RDS file for future use, and launch a local
#' server to explore the results with a browser.
#'
#' @param fpout FastProject object
#' @param ofile the path to save the object in. If NULL, the object is saved
#' in the working directory [default:NULL]
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @return the path of the saved file
#' @aliases saveAndViewResults
#' @export
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
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
#' ## analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- analyze(fp)
#' saveAndViewResults(fp.out)
#' }
setMethod("saveAndViewResults", signature(fpout="FastProject"),
          function(fpout, ofile=NULL, port=NULL, host=NULL, browser=TRUE, name=NULL) {
            if(is.null(ofile)) {
              i <- 1
              ofile <- paste0("./fpout", i, ".rds")
              while (file.exists(ofile)) {
                i <- i+1
                ofile <- paste0("./fpout", i, ".rds")
              }
            }

            saveRDS(fpout, file=ofile)
            viewResults(fpout, port, host, browser, name)
            return(ofile)
          })

#' View results of analysis
#'
#' launch a local server to explore the results with a browser.
#'
#' @param object FastProject object or path to a file containing such an
#' object (saved using saveAndViewResults, or directly using saveRDS)
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @param browser Whether or not to launch the browser automatically (default=TRUE)
#' @param name Name for the sample - is shown at the top of the output report
#' @aliases viewResults
#' @return None
#' @export
#' @rdname viewResults
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
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
#' ## analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- analyze(fp)
#' viewResults(fp.out)
#' }
setMethod("viewResults", signature(object="FastProject"),
          function(object, port=NULL, host=NULL, browser=TRUE, name=NULL) {

            if (!is.null(name)) {
                object@name <- name
            }

            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            launchServer(object, port, host, browser)
          })

#' @rdname viewResults
#' @export
setMethod("viewResults", signature(object="character"),
          function(object, port=NULL, host=NULL, browser=TRUE, name=NULL) {
            fpo <- readRDS(object)
            if(!is(fpo, "FastProject")){
              stop("loaded object not a valid FastProject object")
            }
            viewResults(fpo, port, host, browser, name)
          })

#' create new FastProject object from a subset of the data in an existing one
#' @param fp the FastProject object to subset
#' @param subset the indices of the samples to keep
#' @return a new FastProject object with the new data and the same analysis
#' parameters
createNewFP <- function(fp, subset) {
    .Object <- new("FastProject")
    nexpr <- fp@allData[,subset]
    rownames(nexpr) <- toupper(rownames(nexpr))
    .Object@allData = nexpr
    .Object@exprData <- ExpressionData(nexpr)

    .Object@housekeepingData <- fp@housekeepingData
    .Object@sigData <- fp@sigData

    .Object@precomputedData <- lapply(fp@precomputedData, function(sigscore) {
        sigscore@scores <- sigscore@scores[subset]
        sigscore@sample_labels <- sigscore@sample_labels[subset]
        return(sigscore)
    })

    .Object@weights <- fp@weights[,subset]
    .Object@projection_genes <- fp@projection_genes
    .Object@threshold <- fp@threshold
    .Object@sig_norm_method <- fp@sig_norm_method
    .Object@sig_score_method <- fp@sig_score_method
    .Object@lean = fp@lean
    .Object@perm_wPCA = fp@perm_wPCA
    .Object@pool = fp@pool
    .Object@cellsPerPartition = fp@cellsPerPartition

    return(.Object)
}
