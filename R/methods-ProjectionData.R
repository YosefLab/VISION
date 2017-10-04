#' Wrapper for storing all relevant information for a given projection.
#'
#' Stores a list of Projection objects, filter name, and a logical value indicating whether or not
#' PCA was performed. Also stores clusters, a signature to projection matrix, and relevant gene names
#' and signature / projection keys.

#' Initializes a ProjectionData object for neatly storing all relevant data to a given projection section
#'
#' @param projections List of Projection objects to be stored
#' @param keys Sample names of expression data
#' @param sigProjMatrix Matrix storing the median consistency score per projection, signature pair
#' @param pMatrix Matrix storing the p values for each projection, signature pair
#' @return ProjectionData object
setMethod("initialize", signature(.Object="ProjectionData"),
          function(.Object, projections=NULL, keys, sigProjMatrix, pMatrix, sigClusters) {

            .Object@projections <- projections
            .Object@keys <- keys
            .Object@sigProjMatrix <- sigProjMatrix
            .Object@pMatrix <- pMatrix
            .Object@sigClusters <- sigClusters
            return(.Object)
          }
)

