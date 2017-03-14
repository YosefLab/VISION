#' A wrapper class for clusters computed during the projections step.
#' 
#' Maps a projection name to a list of Cluster Obejcts computed with the projection data.

setMethod("initialize", signature(.Object="ClusterData"), 
          function(.Object, projectionName, clusters) {
            
            .Object@projectionName = projectionName
            .Object@clusters = clusters
            
            return(.Object)
          })