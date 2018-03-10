#' Wrapper for storing all relevant information for a given projection.
#'
#' Stores a list of Projection objects.
#' Also stores clusters, a signature to projection matrix

#' Initializes a ProjectionData object for neatly storing all relevant data to a given projection section
#' Should not be called directly, instead use the `new` syntax
#'
#' @param projections List of Projection objects to be stored
#' @param sigProjMatrix Matrix storing the median consistency score per projection, signature pair
#' @param pMatrix Matrix storing the p values for each projection, signature pair
#' @param sigClusters a list of signature clusters
#' @param treeScore a significance score for the fitted tree
#' @return ProjectionData object
TreeProjectionData <- function(projections=NULL, sigProjMatrix,
                               pMatrix, sigClusters, treeScore) {
    .Object <- new("TreeProjectionData", projections=projections,
                   sigProjMatrix=sigProjMatrix, pMatrix=pMatrix,
                   sigClusters=sigClusters, treeScore=treeScore)
    return(.Object)
    }

