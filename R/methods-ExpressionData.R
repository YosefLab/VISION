
#' Initializes an ExpressionData object
#' 
#' @param data expression data matrix
#' @return ExpressionData object
#' @examples 
#' expr <- readExprtoMatrix("data/expression_matrix.txt")
#' edata <- ExpressionMatrix(expr)
setMethod("initialize", signature(.Object="ExpressionData"),
          function(.Object, data, ...) {
            tryCatch({
              if (missing(data)) {
                data <- matrix(labelDescription = rep(NA, ncol(data)))
              }
            }, error=function(err) {
              stop(conditionMessage(err),
                   "\n ExpressionData 'initialize' could not add data:",
                   "\n perhap data was not entered correctly")
            })
            
            .Object@data <- data
            .Object@fanoFilter <- matrix(NA)
            .Object@thresholdFilter <-matrix(NA)
            .Object@noVarFilter <- matrix(NA)
            return(.Object)
            
          })

#' Prints out expression data
#' 
#' @param object ExpressionData object
#' @return Nothing.
#' @examples 
#' expr <- readExprtoMatrix("data/expression_matrix.txt")
#' edata <- ExpressionMatrix(expr)
#' readExprData(eData) 
setMethod("readExprData", signature("ExpressionData"), function(object) {
	# Prints out the expression data stored in this object.

  print(object@data)
  return()
})

#' Extracts the expression data from the ExpressionData object
#' 
#' @param object ExpressionData object
#' @return Expression data matrix
#' @examples
#' expr <- readExprtoMatrix("data/expression_matrix.txt")
#' edata <- ExpressionMatrix(expr)
#' exprData <- getExprData(eData)
setMethod("getExprData", signature("ExpressionData"), function(object) {
  # Returns the expression data stored in the this object.

  return(object@data)
})

#' Updates the expression data stored in this object
#' 
#' @param object ExpressionData object
#' @param newData new expression data matrix
#' @return ExpressionData object with updated data
setMethod("updateExprData", signature("ExpressionData"), function(object, newData) {
  
  object@data <- newData
  return(object)
})

#' Calculates the specified normalized data matrix
#' 
#'  @param object ExpressionData object
#'  @param func normalization method to apply
#'  @return Normalized data matrix according to function specified.  
#'  @example 
#'  eData <- ExpressionData(expr)
#'  ne <- getNormalizedCopy(eData, "znorm_rows")
setMethod("getNormalizedCopy", signature("ExpressionData"), function(object, func) {

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


