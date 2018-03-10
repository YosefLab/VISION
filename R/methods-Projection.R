## Projection wrapper class

#' Initialize a Projection object
#' Should not be called directly, instead use the `new` syntax
#'
#' @param name Name of the projection
#' @param pData Coordinates of each sample in the projection
#' (NUM_SAMPLES x NUM_COMPONENTS)
#' @param weights a matrix of weights indicatng distances from each point
#' to its closest neighbors
#' @return Projection object

Projection <- function(name, pData=NULL, weights=matrix(NA, 1,1)) {
            .Object <- new("Projection", name=name, pData=pData,
                           weights=weights)
            return(.Object)
            }


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
#' @importFrom stats kmeans
#' @param object Projection object
#' @param method Method by which to cluster the data
#' @param param Parameters for clustering method
#' @return a Cluster object
setMethod("cluster", signature(object = "Projection"),
            function(object, method, param) {

            if(method == "KMeans")
            {
                km <- kmeans(t(object@pData), centers=param)
                clust <- Cluster(method, param, km$centers,
                                    t(as.matrix(km$cluster)))
            }
            else
            {
                stop(paste("Unknown Cluter Method:", method));
            }

            return(clust)
    }
)

#' compute for each vector the weights to apply to it's K nearest neighbors
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom matrixStats rowMaxs
#' @param object the Projecton object
#' @param K number of neughbors to compute ths for
#' @return a weights matrix
setMethod("computeKNNWeights", signature(object = "Projection"),
    function(object, K = 30) {

        if (!is.na(object@weights[1, 1])) {
            return(object@weights)
        }

        n_workers <- getWorkerCount()

        k <- ball_tree_knn(t(object@pData), K, n_workers)
        nn <- k[[1]]
        d <- k[[2]]

        sigma <- rowMaxs(d)
        sparse_weights <- exp(-1 * (d * d) / sigma ^ 2)

        # Normalize row sums = 1
        weightsNormFactor <- rowSums(sparse_weights)
        weightsNormFactor[weightsNormFactor == 0] <- 1.0
        sparse_weights <- sparse_weights / weightsNormFactor

        # load into a sparse matrix
        tnn <- t(nn)
        j <- as.numeric(tnn)
        i <- as.numeric(col(tnn))
        vals <- as.numeric(t(sparse_weights))
        dims <- c(nrow(nn), nrow(nn))
        dimnames <- list(colnames(object@pData), colnames(object@pData))

        weights <- sparseMatrix(i = i, j = j, x = vals,
                                dims = dims, dimnames = dimnames
                               )

        return(weights)
    }
)
