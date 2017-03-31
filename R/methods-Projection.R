## Projection wrapper class


setMethod("initialize", signature(.Object = "Projection"), 
          function(.Object, name, pData=NULL) {
          
            .Object@name = name
            .Object@pData = pData
            
            return(.Object)
          }
)

setMethod("updateProjection", signature(object = "Projection"),
          function(object, data) {
            object@pData <- data
            return(object)
          }
)