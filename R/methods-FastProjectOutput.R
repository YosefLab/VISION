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
#' \dontrun{
#' fpout <- Analyze(fp)
#' saveFPOutAndViewResults(fpout)
#' }
setMethod("saveFPOutAndViewResults", signature(fpout="FastProjectOutput"),
          function(fpout) {
            i <- 1
            ofile <- paste0("FastProject_Output/fpout", i, ".rds")
            while (file.exists(ofile)) {
              i <- i+1
              ofile <- paste0("FastProject_Output/fpout", i, ".rds")
            }

            saveRDS(fpout, file=ofile)
            arg1 <<- ofile

            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            source("FastProject_Output/server_script.R")
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
#' \dontrun{
#' fpout <- Analyze(fp)
#' viewResults(fpout)
#' }
setMethod("viewResults", signature(object="FastProjectOutput"),
          function(object) {

            arg1 <<- object
            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            source("FastProject_Output/server_script.R")
          })
