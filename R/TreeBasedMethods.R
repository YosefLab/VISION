require("igraph")


#' Project the given dataoints onto the tree defined by the vertices (V.pos) and binary adjacency matrix (adj.mat)
#'
#' @param data.pnts (D x N numeric) the spatial coordinates of data points to project
#' @param V.pos (D x K numeric) the spatial coordinates of the tree vertices
#' @param adj.mat (K x K logical) a symmetric binary adjacency matrix (K x K)
#'
#' @details data points are projected on their nearest edge, defined to be the edge that is connected to the nearest node
#' and has minimal length of orthogonal projection. Data points are projected with truncated orthogonal projection:
#' point that fall beyond the edge are projected to the closer node.
#'
#' @return a list of the following:
#' {
#' \item{"spatial"}{The D-dimensional position of the projected data points}
#' \itam{"edge"}{a Nx2 matrix, were line i has the indices identifying the edge that datapoint i was projected on,
#' represented as (node a, node b). For consistency and convenience, it is maintained that a < b}
#' \item{"edge.pos"}{an N-length numeric with values in [0,1], the relative position on the edge of the datapoint.
#' 0 is node a, 1 is node b, .5 is the exact middle of the edge, etc.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(200), nrow = 2)
#' tree <- applySimplePPT(X)
#' proj <- projectOnTree(X, tree[[1]], tree[[2]])
projectOnTree <- function(data.pnts, V.pos, adj.mat) {
  # find closest principle point
  distmat <- sqdist(t(data.pnts), t(V.pos))
  major.bool <- t(apply(distmat, 1, function(x) {x == min(x)}))
  major.ind <- apply(major.bool, 1, which)

  # find all edges connected to closest principle point
  distmat[major.bool] <- NA # replace closest with NA
  neigh <- adj.mat[major.ind,] # get neighbors of nearest pp
  distmat[neigh == 0] <- NA # remove non-neighbors
  projections <- sapply(1:NCOL(data.pnts), function(i) {
    # for each datapoint, find edge with smallest orthogonal projection
    edges <- t(apply(cbind(major.ind[i], which(!is.na(distmat[i,]))), 1, sort))

    if(NROW(edges) > 1) { ## Not a leaf
      edge.p1 <- V.pos[,edges[,1]]
      edge.p2 <- V.pos[,edges[,2]]

      line <- edge.p2 - edge.p1

      pos <- colSums((data.pnts[,i] - edge.p1) * line) / colSums(line ^ 2)
      rpos <- pmax(0, pmin(1, pos)) ## relative positin on the edge
      spos <- edge.p1 + t(rpos * t(line)) ## spatial position of projected points

      # the best edge is the one with the shortest projection
      best <- which.min(sqrt(colSums((data.pnts[,i] - spos) ^ 2)))

      return(c(edges[best,], rpos[best], spos[,best]))
    } else { # closest node is a leaf, only one appropriate edge

      edge.p1 <- V.pos[,edges[1]]
      edge.p2 <- V.pos[,edges[2]]

      line <- edge.p2 - edge.p1

      pos <- sum((data.pnts[,i] - edge.p1) * line) / sum(line ^ 2)
      rpos <- pmax(0, pmin(1, pos))
      spos <- edge.p1 + t(rpos * t(line))

      return(c(edges, rpos, spos))
    }
  })

  return(list(edges = projections[1:2,],
              edges.pos = projections[3,],
              spatial = projections[-c(1:3),]))
}

calculateTreeKNN <- function(princ.pnts, adj.mat, proj.edges, proj.edge.pos, K) {
  # create principle graph
  # create full graph (with projected data points)
  # Traverse with BFS on the graph and compute distances from each node to next node

}
