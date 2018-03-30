#' Wrapper for storing all relevant information for a given projection.
#'
#' Stores a list of Projection objects
#' Also stores clusters, a signature to projection matrix

#' Initializes a ProjectionData object for neatly storing all relevant data to a given projection section
#'
#' @param projections List of Projection objects to be stored
#' @param sigProjMatrix Matrix storing the median consistency score per projection, signature pair
#' @param pMatrix Matrix storing the p values for each projection, signature pair
#' @param sigClusters the signature clusters
#' @return ProjectionData object
ProjectionData <- function(sigProjMatrix, pMatrix, sigClusters, emp_pMatrix) {

            .Object <- new("ProjectionData", sigProjMatrix = sigProjMatrix,
                           pMatrix = pMatrix, sigClusters = sigClusters,
                           emp_pMatrix = emp_pMatrix)

            return(.Object)
            }
