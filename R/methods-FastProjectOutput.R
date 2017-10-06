#' Initialize a new FastProjectOuptput object
#' should not be called directly, use the `new` syntax
#'
#' @param .Object an object
#' @param eData ExpressionData matrix
#' @param filterModuleList List of FilterModuleData objects
#' @param sigMatrix Matrix of signature scores (NUM_SIGS x NUM_CELLS)
#' @param sigList List of Signatures
#' @param fpParams List of the FastProject parameters that were used to
#' generate this object
#' @param pools a list of pooling information for pooling single cells into
#' representative super-cell clusters
#' @return FastProjectOutput object.
setMethod("initialize", signature(.Object="FastProjectOutput"),
            function(.Object, eData, filterModuleList, sigMatrix, sigList,
                    fpParams, pools) {

            .Object@exprData <- eData
            .Object@filterModuleList <- filterModuleList
            .Object@sigMatrix <- sigMatrix
            .Object@sigList <- sigList
            .Object@fpParams <- fpParams
            .Object@pools <- pools

            return(.Object)
            }
)

#' Save the FastProjectOutput object as an .RDS file and view the results on a
#' localhost
#'
#' Save the results object as an RDS file for future use, and launch a local
#' server to explore the results with a browser.
#'
#' @param fpout FastProjectOutput object
#' @return None
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
#' fp.out <- Analyze(fp)
#' saveFPOutAndViewResults(fp.out)
#' }
setMethod("saveFPOutAndViewResults", signature(fpout="FastProjectOutput"),
            function(fpout) {
            i <- 1
            ofile <- paste0("./fpout", i, ".rds")
            while (file.exists(ofile)) {
                i <- i+1
                ofile <- paste0("./fpout", i, ".rds")
            }

            saveRDS(fpout, file=ofile)
            arg1 <<- ofile

            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            path <- find.package("FastProjectR")
            source(paste0(path, "/FastProjectR_Output/server_script.R"))
            })

#' View results of analysis without saving output object
#'
#' launch a local server to explore the results with a browser.
#'
#' @param object FastProjectOutput object
#' @aliases viewResults
#' @return None
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
#'                      housekeeping = hkg,
#'                      signatures = sigs)
#'
#' ## Analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- Analyze(fp)
#' viewResults(fp.out)
#' }
setMethod("viewResults", signature(object="FastProjectOutput"),
            function(object) {

            arg1 <<- object
            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            path <- find.package("FastProjectR")
            source(paste0(path, "/FastProjectR_Output/server_script.R"))
            })
