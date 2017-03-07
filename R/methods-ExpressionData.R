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

setMethod("readExprData", signature("ExpressionData"), function(object) {
  print(object@data)
  return()
})

setMethod("getExprData", signature("ExpressionData"), function(object) {
  return(object@data)
})

setMethod("updateExprData", signature("ExpressionData"), function(object, newData=data.frame()) {
  object@data <- newData
  return(object)
})

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


