
#' Create a TreeProjection object
#' Should not be called directly, instead use the `new` syntax
#'
#' @param pData the projected coordinates of the data
#' @param name the name of the projection
#' @param vData the projected coordinates of the tree nodes
#' @param adjMat the adjacency matrix of the tree
#' @return a TreeProjection object
TreeProjection <- function(pData, name, vData, adjMat) {
    proj <- projectOnTree(data.pnts = pData,
                            V.pos = vData,
                            princAdj = adjMat)
    .Object <- new("TreeProjection", pData=pData, name=name,
                   vData=vData, adjMat=adjMat, edgeAssoc=proj$edges,
                   edgePos=proj$edgePos)
    return(.Object)
    }


#' Compute KNN weights based on geodesic distances for TreeProjection objects
#'
#' @param object a TreeProjection object
#' @param K the number of nearest neighbors to look at
#' @param BPPARAM the parallelization backen to use
#' @return an all-pars distance matrix
setMethod("computeKNNWeights", signature(object = "TreeProjection"),
            function(object, K=30, BPPARAM=bpparam()) {
            distmat <- calculateTreeDistances(princPnts = object@vData,
                                                princAdj = object@adjMat,
                                                edgeAssoc = object@edgeAssoc,
                                                edgePos = object@edgePos)

            kQuantile <- K / NCOL(object@pData)
            knnmat <- apply(distmat, 1, function(d) {
                partition <- stats::quantile(d, kQuantile)
                d[d > partition] <- Inf
                return(d)
            })

            sigma <- apply(knnmat, 1, max)
            weights <- exp(-1 * (knnmat * knnmat) / sigma^2)

            weights[is.na(weights)] <- 0.0

            weightsNormFactor <- Matrix::rowSums(weights)
            weightsNormFactor[weightsNormFactor == 0] <- 1.0
            weightsNormFactor[is.na(weightsNormFactor)] <- 1.0
            weights <- weights / weightsNormFactor

            return(weights)

            })
