#' Compute KNN weights based on geodesic distances for Trajectory objects
#' @importFrom stats quantile
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom matrixStats rowMaxs
#' @param object a Trajectory object
#' @param K the number of nearest neighbors to look at
#' @return a list of two items:
#'          indices: matrix, cells X neighbors
#'              Each row specifies indices of nearest neighbors
#'          weights: matrix, cells X neighbors
#'              Corresponding weights to nearest neighbors
setMethod("computeKNNWeights", signature(object = "Trajectory"),
            function(object, K = round(sqrt(nrow(object@pData))) ) {

            # Todo: Fill in stuff here
            # Similar to the one in methods-TreeProjection

            return(list(indices = nn, weights = sparse_weights))
            })
