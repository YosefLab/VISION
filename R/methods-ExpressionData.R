setMethod("initialize", signature(.Object="ExpressionData"),
          function(.Object, data=data.frame(), ...) {
            tryCatch({
              if (missing(data)) {
                data <- data.frame(labelDescription = rep(NA, ncol(data)))
              }
            }, error=function(err) {
              stop(conditionMessage(err),
                   "\n ExpressionData 'initialize' could not add data:",
                   "\n perhap data was not entered correctly")
            })
            return(.Object)
            
          })

setMethod("readExprData", signature("ExpressionData"), function(object) {
  print(object@data)
  return()
})

setMethod("updateExprData", signature("ExpressionData"), function(object, newData=data.frame()) {
  object@data <- newData
  return()
})


