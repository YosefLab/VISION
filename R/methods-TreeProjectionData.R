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
#' @param emp_pMatrix Matrix storing the empirical p values for each projection, signature pair
#' @param sigClusters a list of signature clusters
#' @param treeScore a significance score for the fitted tree
#' @return ProjectionData object
TreeProjectionData <- function(latentTree, projections, sigProjMatrix,
                               pMatrix, sigClusters, treeScore, emp_pMatrix) {

    .Object <- new("TreeProjectionData", latentTree = latentTree,
                   projections = projections, sigProjMatrix = sigProjMatrix,
                   pMatrix = pMatrix, sigClusters = sigClusters,
                   treeScore = treeScore, emp_pMatrix = emp_pMatrix)

    return(.Object)
    }

