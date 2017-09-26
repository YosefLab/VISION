## Projection wrapper class
require("cluster")

#' Initialize a Projection object
#'
#' @param name Name of the projection
#' @param pData Coordinates of each sample in the projection (NUM_SAMPLES x NUM_COMPONENTS)
#' @return Projection object
setMethod("initialize", signature(.Object = "Projection"),
          function(.Object, name, pData=NULL, ppt_c=matrix(NA, 1, 1)) {

            .Object@name = name
            .Object@pData = pData
            .Object@PPT_C = ppt_c

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


setMethod("computeKNNWeights", signature(object = "Projection"),
          function(object, K = 30, numCores = 1) {
            weights <- matrix(0L, nrow=NCOL(proj@pData), ncol=NCOL(proj@pData))
            k <- ball_tree_knn(t(proj@pData), K, numCores)
            nn <- k[[1]]
            d <- k[[2]]

            sigma <- apply(d, 1, max)
            sparse_weights <- exp(-1 * (d * d) / sigma^2)

            weights <- load_in_knn(nn, sparse_weights)
            #d <- dist.matrix(t(proj@pData), method="euclidean")
            #weights <- exp( (-1 * (d * d)) / (NEIGHBORHOOD_SIZE)^2)

            weightsNormFactor <- Matrix::rowSums(weights)
            weightsNormFactor[weightsNormFactor == 0] <- 1.0
            weightsNormFactor[is.na(weightsNormFactor)] <- 1.0
            weights <- weights / weightsNormFactor

            return(weights)
          }
)

