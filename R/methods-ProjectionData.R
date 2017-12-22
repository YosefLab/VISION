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
#' @param sigClusters the signature clusters
#' @return ProjectionData object
ProjectionData <- function(projections=NULL, keys, sigProjMatrix,
                           pMatrix, sigClusters, emp_pMatrix) {
            .Object <- new("ProjectionData", projections=projections,
                           keys=keys, sigProjMatrix=sigProjMatrix,
                           pMatrix=pMatrix, sigClusters=sigClusters, emp_pMatrix=emp_pMatrix)
            return(.Object)
            }

