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

  # minors <- row.argmin(distmat)
  minor <- apply(distmat, 1, which.min)

  edges <- t(apply(cbind(major.ind, minor), 1, sort))

  # get spatial positions of the edges
  edge.p1 <- V.pos[,edges[,1]]
  edge.p2 <- V.pos[,edges[,2]]

  line <- edge.p2 - edge.p1

  # relative position on the edge
  score <- pmax(0, pmin(1, colSums((data.pnts - edge.p1) * line) / colSums(line ^ 2)))

  #spatial position of the projected point
  pos <- edge.p1 + t(score * t(line))

  return(list("spatial" = pos, "edge" = edges, "edge.pos" = score))
}

#' Get pseudotime alignment of cells on the path from v.s to v.t
#'
#' @param adj.mat
getPseudoTime <- function(adj.mat, edges, edges.pos, v.s, v.t) {

}
