#' Initialize a new TrajectoryProjection object.
#'
#' This represents a 2d view of data on a projection
#'
#' @param name label for this projection
#' @param vData #Milestones x 2 matrix with coordinates for milestones
#' @param pData #Cells x 2 matrix with coordinates for cells
#' @param adjMat #Milestones x #Milestones adjacency matrix (0/1)
#' @return TrajectoryProjection object
TrajectoryProjection <- function(name, vData, pData, adjMat) {

    if (!ncol(vData) == 2){
        stop("bad vData (not 2 columns)")
    }
    if (!ncol(pData) == 2){
        stop("bad pData (not 2 columns)")
    }
    if (nrow(adjMat) != ncol(adjMat)){
        stop("adjMat should be symmetric")
    }
    if (nrow(adjMat) != nrow(vData)){
        stop("nrow(vData) should = nrow(adjMat)")
    }

    .Object <- new("TrajectoryProjection", name = name,
                   vData = vData, pData = pData, adjMat = adjMat)

    return(.Object)
}
