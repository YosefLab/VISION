
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
#' @importFrom stats quantile
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom matrixStats rowMaxs
#' @param object a TreeProjection object
#' @param K the number of nearest neighbors to look at
#' @return an all-pars distance matrix
setMethod("computeKNNWeights", signature(object = "TreeProjection"),
            function(object, K=30) {
            distmat <- calculateTreeDistances(princPnts = object@vData,
                                                princAdj = object@adjMat,
                                                edgeAssoc = object@edgeAssoc,
                                                edgePos = object@edgePos)

            kQuantile <- K / NCOL(object@pData)
            knnmat <- apply(distmat, 1, function(d) {
                partition <- quantile(d, kQuantile)
                d[d > partition] <- Inf
                return(d)
            })

            nn <- t(apply(distmat, 1, function(r) {
                            order(r)[1:K]
            }))

            d <- lapply(seq(nrow(nn)), function(i) {
                            distmat[i, nn[i, ]]
            })
            d <- do.call(rbind, d)

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

            weights <- sparseMatrix(i = i, j = j, x = vals,
                                    dims = c(nrow(nn), nrow(nn))
                                    )

            return(weights)

            })
