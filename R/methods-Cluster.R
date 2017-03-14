#' Wrapper class for a particular cluster. 
#' 
#' Maps a cluster type to the the resulting cluster data.

setMethod("initialize", signature(.Object="Cluster"),
          function(.Object, name, centers, data) {
            
            .Object@name <- name
            .Object@centers <- centers
            .Object@data <- data
            
            return(.Object)
          })