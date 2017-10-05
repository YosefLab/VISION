#' Wrapper class for a particular cluster.
#' Should be called using the `new` syntax.
#' Maps a cluster type to the the resulting cluster data.
#' @param .Object An object.
#' @param method clustering method used
#' @param param clusterng parameter used
#' @param centers cluster centroid
#' @param data cluster-assocated data
#' @return a Cluster object
setMethod("initialize", signature(.Object="Cluster"),
          function(.Object, method, param, centers, data) {

            .Object@method <- method
            .Object@param <- param
            .Object@centers <- centers
            .Object@data <- data

            return(.Object)
          })
