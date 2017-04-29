## Projection wrapper class
require("cluster")


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

setMethod("cluster", signature(object = "Projection"),
          function(object, method, param) {

            if(method == "KMeans")
            {
                km <- kmeans(t(object@pData), centers=param)
                clust <- Cluster(method, param, km$centers, t(as.matrix(km$cluster)))
            }
            else
            {
                stop(paste("Unknown Cluter Method:", method));
            }

            return(clust)
    }
)
