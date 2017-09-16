require("igraph")


#' Project the given dataoints onto the tree defined by the vertices (V.pos) and binary adjacency matrix (adj.mat)
#'
#' @param data.pnts the data points to project (D x N)
#' @param V.pos the positions of the tree vertices (D x K)
#' @param adj.mat a binary, symmetric adjacency matrix (K x K)
#'
#' @details the points are projected in two stages. First, the closest edge is found, by finding the closest vertex,
#' and then finding the neighbor of that vertex that is closest to the datapoint. The edge connecting the two vertices
#' is considered the closest edge. Second, the datapoint is linearly projected onto the line of the edge, in a truncated manner,
#' so points that lie beyond to bounds of the edge itself are projected onto the vertex.
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

  # find closest neighbor of the closest principle point
  distmat[major.bool] <- NA # replace closest with NA
  neigh <- adj.mat[major.ind,] # get neighbors of nearest pp
  distmat[neigh == 0] <- NA # remove non-neighbors
  projections <- sapply(1:NCOL(data.pnts), function(i) {
    edges <- t(apply(cbind(major.ind[i], which(!is.na(distmat[i,]))), 1, sort))

    if(NROW(edges) > 1) { # find the best edge to project on, based
      edge.p1 <- V.pos[,edges[,1]]
      edge.p2 <- V.pos[,edges[,2]]

      line <- edge.p2 - edge.p1

      pos <- colSums((data.pnts[,i] - edge.p1) * line) / colSums(line ^ 2)
      rpos <- pmax(0, pmin(1, pos))
      spos <- edge.p1 + t(rpos * t(line))

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
#
# projectOnLine <- function(dp, ep1, ep2) {
#   line <- ep2 - ep1
#
#   # relative position on the edge
#   pos <- 1, sum((dp - ep1) * line) / sum(line ^ 2)
#
#   rpos <- pmax(0, pmin(pos))
#
#   #spatial position of the projected point
#   spos <- ep1 + t(rpos * t(line))
# }
#
#
# getPseudoTime <- function(adj.mat, edges, edges.pos, v.s, v.t) {
#
# }
#
# edge.orders <- function(princ.pnts, edges, edges.pos) {
#
# }
#
# edge.order <- function(edges, edges.pos, from, to, len.factor = 1) {
#   key <- sort(c(from, to))
#
# }
