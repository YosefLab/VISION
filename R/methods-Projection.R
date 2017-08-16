## Projection wrapper class
require("cluster")

#' Initialize a Projection object
#' 
#' @param name Name of the projection
#' @param pData Coordinates of each sample in the projection (NUM_SAMPLES x NUM_COMPONENTS)
#' @return Projection object
setMethod("initialize", signature(.Object = "Projection"), 
          function(.Object, name, pData=NULL) {
          
            .Object@name = name
            .Object@pData = pData
            
            return(.Object)
          }
)

#' Updates the coordinate data stored in this object.
#' 
#' @param object Projection object
#' @param data New data to be stored in the object
#' @return Updated Projection object. 
setMethod("updateProjection", signature(object = "Projection"),
          function(object, data) {
            object@pData <- data
            return(object)
          }
)

#' Clusters the projection according to some method
#' 
#' @param object Projection object
#' @param method Method by which to cluster the data
#' @param param Parameters for clustering method
#' @examples 
#' p <- Projection("PCA", pData)
#' cl <- cluster(p, "KMeans", 10)
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
