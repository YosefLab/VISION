#' Wrapper class for a particular cluster.
#' Maps a cluster type to the the resulting cluster data.
#' @param method clustering method used
#' @param param clusterng parameter used
#' @param centers cluster centroid
#' @param data cluster-assocated data
#' @return a Cluster object
Cluster <- function(method, param, centers, data) {
            .Object <- new("Cluster", method=method, param=param,
                           centers=centers, data=data)
            return(.Object)
            }
