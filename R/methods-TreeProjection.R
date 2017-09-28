

setMethod("initialize", signature(.Object="TreeProjection"),
          function(.Object, ..., vData, adjMat) {
            .Object@vData = vData
            .Object@adjMat = adjMat

            callNextMethod(.Object, ...)

            proj <- projectOnTree(data.pnts = .Object@pData,
                                  v.pos = object.vData,
                                  princAdj = object.adjMat)
            .Object@edgeAssoc <- proj$edges
            .Object@edgePos <- proj$edgePos
          }
          )


#' Compute KNN weights based on geodesic distances for TreeProjection objects
#'
#' @param object a TreeProjection object
#' @param K the number of nearest neighbors to look at
#' @param numCores the number of cores to utilize
#'
setMethod("computeKNNWeights", signature(object = "TreeProjection"),
          function(object, K=30, numCores=1) {
            distmat <- calculateTreeDistances(princPnts = object@vData,
                                              princAdj = object@adjMat,
                                              edgeAssoc = proj$edges,
                                              edgePos = proj$edgesPos)

            kQuantile <- K / NCOL(object@pData)
            knnmat <- apply(distmat, 1, function(d) {
              partition <- quantile(d, kQuantile)
              d[d > partition] <- Inf
              return(d)
            })

            sigma <- apply(knnmat, 1, max)
            weights <- exp(-1 * (knnmat * knnmat) / sigma^2)

            weightsNormFactor <- Matrix::rowSums(weights)
            weightsNormFactor[weightsNormFactor == 0] <- 1.0
            weightsNormFactor[is.na(weightsNormFactor)] <- 1.0
            weights <- weights / weightsNormFactor

          })
