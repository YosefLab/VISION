#' Wrapper class for a particular cluster. 
#' 
#' Maps a cluster type to the the resulting cluster data.

setMethod("initialize", signature(.Object="Cluster"),
          function(.Object, method, param, centers, data) {
            
            .Object@method <- method
            .Object@param <- param
            .Object@centers <- centers
            .Object@data <- data
            
            return(.Object)
          })
