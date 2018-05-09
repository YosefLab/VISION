#' Initialize a new Trajectory object.
#'
#' @param input trajectory to model cell progression.  Wrapped result
#' of a trajectory inference by the dynverse/dynwrap library
#' @return Trajectory object
Trajectory <- function(input) {

    # Create the adjacency matrix
    network <- as.data.frame(input$milestone_network)
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

    .Object <- new("Trajectory", adjMat = adjMat, progressions = progressions)

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
            function(object, K = round(sqrt(nrow(object@progressions))) ) {

            edgePos <- object@progressions$position
			names(edgePos) <- rownames(object@progressions)
			edgeAssoc <- t(object@progressions[, c("from", "to")])

            distmat <- calculateTrajectoryDistances(adjMat = object@adjMat,
                                                edgeAssoc = edgeAssoc,
                                                edgePos = edgePos)

            kQuantile <- K / nrow(object@progressions)
            knnmat <- apply(distmat, 1, function(d) {
                partition <- quantile(d, kQuantile, na.rm=T)
                d[d > partition] <- Inf
                return(d)
            })

            nn <- t(apply(distmat, 1, function(r) {
                            order(r)[1:K]
            }))

            d <- lapply(seq(nrow(nn)), function(i) {
                            distmat[i, nn[i, ]]
            })
            d <- do.call(rbind, d)

	    d[is.na(d)] = 0
            sigma <- rowMaxs(d)
            sigma[sigma == 0] <- 1.0 # occurs if all nearest neighbors at same point
            sparse_weights <- exp(-1 * (d * d) / sigma ^ 2)

            # Normalize row sums = 1
            weightsNormFactor <- rowSums(sparse_weights)
            weightsNormFactor[weightsNormFactor == 0] <- 1.0
            sparse_weights <- sparse_weights / weightsNormFactor

            rownames(nn) <- rownames(object)
            rownames(d) <- rownames(object)

            return(list(indices = nn, weights = sparse_weights))
            })


#' Generate meta-data associated with this trajectory
#'
#' Creates a categorical variable mapping cells to edges
#' and numeric variables for their position along edges
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph ends
#' @importFrom igraph E
#'
#' @param trajectory Trajectory on which to operate
#' @return metaData dataframe with meta-data
createTrajectoryMetaData <- function(trajectory){

    adjMat <- trajectory@adjMat
    progressions <- trajectory@progressions
    net <- igraph::graph_from_adjacency_matrix(adjMat, weighted = TRUE,
                                               mode = "undirected")
    edges <- igraph::ends(net, igraph::E(net), names = TRUE)

    meta <- data.frame(row.names = rownames(progressions))
    meta[, "TrajectoryEdge"] <- ""

    for (i in seq(nrow(edges))) {
        from <- edges[i, 1]
        to <- edges[i, 2]

        cells <- progressions[
                     (progressions$from == from) & (progressions$to == to),
                     , drop = FALSE
                     ]

        cells_i <- progressions[
                       (progressions$from == to) & (progressions$to == from),
                       , drop = FALSE
                       ]

        cells_i$position <- 1 - cells_i$position
        cells <- rbind(cells, cells_i)

        edge_name <- paste(from, to, sep = "->")
        position_var <- paste0("Position: ", edge_name)

        meta[rownames(cells), "TrajectoryEdge"] <- edge_name
        meta[, position_var] <- NA
        meta[rownames(cells), position_var] <- cells$position
    }

    meta$TrajectoryEdge <- as.factor(meta$TrajectoryEdge)

    return(meta)

}



#' Generate 2d representations of a trajectory model
#'
#' Creates 2d layouts of the milestone network and translates the
#' cell positions into these layouts
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph ends
#' @importFrom igraph layout_with_fr
#' @importFrom igraph layout_with_dh
#' @importFrom igraph layout_as_tree
#' @importFrom igraph layout_with_mds
#' @importFrom igraph E
#' @importFrom igraph V
#'
#' @param trajectory Trajectory on which to operate
#' @return trajectoryProjections list of TrajectoryProjection
generateTrajectoryProjections <- function(trajectory) {

    progressions <- trajectory@progressions
    adjMat <- trajectory@adjMat

    net <- igraph::graph_from_adjacency_matrix(adjMat, weighted = TRUE,
                                               mode = "undirected")


    adjMatInv <- adjMat ** -1
    adjMatInv[is.infinite(adjMatInv)] <- 0

    # some algorithms use weights to represent 'inverse' distance
    invnet <- igraph::graph_from_adjacency_matrix(adjMatInv, weighted = TRUE,
                                               mode = "undirected")

    edges <- igraph::ends(net, igraph::E(net), names = TRUE)
    adjMatBinary <- (adjMat > 0) * 1

    tp_list <- list()

    # layout with Fruchterman-Reingold algorithm
    vData <- igraph::layout_with_fr(invnet)
    rownames(vData) <- igraph::V(net)$name

    pData <- translateCellPositions(progressions, vData, edges)

    tp <- TrajectoryProjection(name = "FR", pData = pData, vData = vData,
                               adjMat = adjMatBinary)

    tp_list <- c(tp_list, tp)

    # layout with Davidson-Harel
    vData <- igraph::layout_with_dh(net)
    rownames(vData) <- igraph::V(net)$name

    pData <- translateCellPositions(progressions, vData, edges)

    tp <- TrajectoryProjection(name = "DH", pData = pData, vData = vData,
                               adjMat = adjMatBinary)

    tp_list <- c(tp_list, tp)

    # layout with Davidson-Harel
    vData <- igraph::layout_as_tree(net)
    rownames(vData) <- igraph::V(net)$name

    pData <- translateCellPositions(progressions, vData, edges)

    tp <- TrajectoryProjection(name = "Tree", pData = pData, vData = vData,
                               adjMat = adjMatBinary)

    tp_list <- c(tp_list, tp)

    # layout with MDS
    vData <- igraph::layout_with_mds(net)
    rownames(vData) <- igraph::V(net)$name

    pData <- translateCellPositions(progressions, vData, edges)

    tp <- TrajectoryProjection(name = "MDS", pData = pData, vData = vData,
                               adjMat = adjMatBinary)

    tp_list <- c(tp_list, tp)

    names(tp_list) <- vapply(tp_list, function(x) x@name, FUN.VALUE = "")

    return(tp_list)
}


#' Translate cell positions
#'
#' Maps cell positions between edges
#'
#' @param progressions data.frame describing cell positions between milestones
#' @param vData Mx2 matrix mapping miletones into 2d
#' @param edges Edges x 2 matrix describing connectivity between edges
#' @return pData Cells x 2 matrix with cell positions in 2d
translateCellPositions <- function(progressions, vData, edges) {
    pData <- lapply(seq(nrow(edges)), function(i) {
        from <- edges[i, 1]
        to <- edges[i, 2]

        edge_dist <- sum( (vData[from, ] - vData[to, ]) ** 2) ** .5


        cells <- progressions[
                     (progressions$from == from) & (progressions$to == to),
                     , drop = FALSE
                     ]

        cells_i <- progressions[
                       (progressions$from == to) & (progressions$to == from),
                       , drop = FALSE
                       ]

        cells_i$position <- 1 - cells_i$position
        cells <- rbind(cells, cells_i)

        coordinates <- matrix(numeric(nrow(cells) * 2),
                              nrow = nrow(cells), ncol = 2,
                              dimnames = list(rownames(cells), c("x", "y")))

        coordinates[, "x"] <- cells$position * edge_dist
        coordinates[, "y"] <- rnorm(nrow(coordinates), sd = .1)

        dx <- vData[to, 1] - vData[from, 1]
        dy <- vData[to, 2] - vData[from, 2]

        sinTheta <- dy / edge_dist
        cosTheta <- dx / edge_dist

        R <- matrix(c(cosTheta, sinTheta, -1 * sinTheta, cosTheta), nrow = 2)

        coordinates <- coordinates %*% t(R) # rotation

        coordinates <- t(t(coordinates) + vData[from, ]) # offset

        return(coordinates)
    })

    pData <- do.call(rbind, pData)
    pData <- pData[rownames(progressions), , drop = FALSE]
    return(pData)
}
