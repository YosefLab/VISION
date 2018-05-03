#' Initialize a new Trajecotry object.
#'
#' @param input trajectory to model cell progression.  Wrapped result
#' of a trajectory inference by the dynverse/dynwrap library
#' @return Trajectory object
Trajectory <- function(input) {

    # Create the adjacency matrix
    network <- input$milestone_network
    milestone_ids <- union(unique(network$from), unique(network$to))
    adjMat <- matrix(0.0, nrow = length(milestone_ids), ncol = length(milestone_ids),
                     dimnames = list(milestone_ids, milestone_ids))

    for (i in seq(nrow(network))){
        from <- as.character(network[i, "from"])
        to <- as.character(network[i, "to"])
        edgeLength <- network[i, "length"]

        adjMat[from, to] <- edgeLength
        adjMat[to, from] <- edgeLength
    }

    progressions <- input$progressions
    rownames(progressions) <- progressions$cell_id
    progressions$cell_id <- NULL

    colnames(progressions) <- gsub("percentage", "position",
                                   colnames(progressions))

    if (!("from" %in% colnames(progressions))){
        stop("input Trajectory missing 'from' column")
    }

    if (!("to" %in% colnames(progressions))){
        stop("input Trajectory missing 'to' column")
    }

    if (length(setdiff(unique(progressions$from), milestone_ids)) > 0) {
        stop("milestones in progressions$from don't match those in milestone_network")
    }

    if (length(setdiff(unique(progressions$to), milestone_ids)) > 0) {
        stop("milestones in progressions$to don't match those in milestone_network")
    }

    .Object <- new("Signature", adjMat = adjMat, progressions = progressions)

    return(.Object)
}

#' Compute KNN weights based on geodesic distances for Trajectory objects
#' @importFrom stats quantile
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom matrixStats rowMaxs
#' @param object a Trajectory object
#' @param K the number of nearest neighbors to look at
#' @return a list of two items:
#'          indices: matrix, cells X neighbors
#'              Each row specifies indices of nearest neighbors
#'          weights: matrix, cells X neighbors
#'              Corresponding weights to nearest neighbors
setMethod("computeKNNWeights", signature(object = "Trajectory"),
            function(object, K = round(sqrt(nrow(object@pData))) ) {

            # Todo: Fill in stuff here
            # Similar to the one in methods-TreeProjection

            return(list(indices = nn, weights = sparse_weights))
            })
