
#' Initializes an ExpressionData object
#' should not be called directly, use the `new` syntax
#' @param data expression data matrix
#' @param ... additional arguments
#' @return ExpressionData object
ExpressionData <- function(data, ...) {
            .Object <- new("ExpressionData", data=data)
            return(.Object)
            }

#' Extracts the expression data from the ExpressionData object
#'
#' @param object ExpressionData object
#' @return Expression data matrix
setMethod("getExprData", signature("ExpressionData"), function(object) {
    # Returns the expression data stored in the this object.

    return(object@data)
})

#' Calculates the specified normalized data matrix
#' @param object ExpressionData object
#' @param func normalization method to apply
#' @return Normalized data matrix according to function specified.
setMethod("getNormalizedCopy", signature("ExpressionData"),
            function(object, func) {

    if (func == "none") {
    return(noNormalization(object@data))
    } else if (func == "znorm_columns") {
    return(colNormalization(object@data))
    } else if (func == "znorm_rows") {
    return(rowNormalization(object@data))
    } else if (func == "znorm_rows_then_columns") {
    return(rowAndColNormalization(object@data))
    } else if (func == "rank_norm_columns") {
    return(colRankNormalization(object@data))
    }
    stop("Normalization method not recognized.")
})


