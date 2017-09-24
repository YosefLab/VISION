require("igraph")


#' Project the given dataoints onto the tree defined by the vertices (V.pos) and binary adjacency matrix (princAdj)
#'
#' @param data.pnts (D x N numeric) the spatial coordinates of data points to project
#' @param V.pos (D x K numeric) the spatial coordinates of the tree vertices
#' @param princAdj (K x K logical) a symmetric binary adjacency matrix (K x K)
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
projectOnTree <- function(data.pnts, V.pos, princAdj) {
  # find closest principle point
  distmat <- sqdist(t(data.pnts), t(V.pos))
  major.bool <- t(apply(distmat, 1, function(x) {x == min(x)}))
  major.ind <- apply(major.bool, 1, which)

  # find all edges connected to closest principle point
  distmat[major.bool] <- NA # replace closest with NA
  neigh <- princAdj[major.ind,] # get neighbors of nearest pp
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

calculateTreeDistances <- function(princPnts, princAdj, edgeAssoc, edgePos) {
  # get all distances in principle tree
  princAdjW <- sqdist(t(princPnts), t(princPnts)) * princAdj

  princGraph <- graph_from_adjacency_matrix(princAdjW, weighted = T, mode = "undirected")
  nodeDistmat <- distances(princGraph)

  princEdges <- apply(get.edgelist(princGraph), 1, as.numeric)
  edgeToPnts <- apply(princEdges, 2, function(x) { apply(edgeAssoc==x, 2, all) })

  distmat <- matrix(rep(NA, NROW(edgeToPnts) ^ 2), NROW(edgeToPnts))

  # compute intra-edge distances. Store in list for later calclations and set values in result matrix
  intraDist <- list()
  for (i in 1:NCOL(princEdges)) {
    inEdgeDist <- calcIntraEdgeDistMat(edge.len = nodeDistmat[princEdges[1,i],
                                                              princEdges[2,i]],
                                       edge.pos = edgePos[edgeToPnts[,i]])
    intraDist[[i]] <- inEdgeDist[,c(1,NCOL(inEdgeDist))]
    distmat[edgeToPnts[,i], edgeToPnts[,i]] <- inEdgeDist[-c(1,NROW(inEdgeDist)), -c(1,NROW(inEdgeDist))]
  }

  # for each pair of edges, calculate inter-edge distances and set them in the empty matrix
  for (i in 1:(NCOL(princEdges)-1)) {
    for (j in (i+1):NCOL(princEdges)) {
      ## figure out which pair is the right one (one with shortest distance)
      edge1NodeInd <- which.min(nodeDistmat[princEdges[,i],princEdges[2,j]])
      edge2NodeInd <- which.min(nodeDistmat[princEdges[edge1NodeInd, i], princEdges[,j]])
      pathLength <- nodeDistmat[princEdges[edge1NodeInd,i], princEdges[edge2NodeInd,j]]

      ## get corresponding distance vectors from intraList
      edge1DistMat <- intraDist[[i]]
      edge1DistVec <- edge1DistMat[-c(1,NROW(edge1DistMat)), edge1NodeInd]
      edge2DistMat <- intraDist[[j]]
      edge2DistVec <- edge2DistMat[-c(1,NROW(edge2DistMat)), edge2NodeInd]

      ## set interedge distances in matrix
      interDistmat <- calcInterEdgeDistMat(v1.dist = edge1DistVec,
                                           v2.dist = edge2DistVec,
                                           path.length = pathLength)
      distmat[edgeToPnts[,i], edgeToPnts[,j]] <- interDistmat
    }
  }

  # since we don't set all coordinates, but the matrix is symmetric
  return(pmax(distmat, t(distmat), na.rm = T))
}

calcIntraEdgeDistMat <- function(edge.len, edge.pos) {
  edge.pos <- c(0, edge.pos, 1) * edge.len
  pos.mat <- replicate(length(edge.pos), edge.pos)
  return(abs(pos.mat - t(pos.mat)))
}

calcInterEdgeDistMat <- function(v1.dist, v2.dist, path.length) {
  # note that input vector must not contain distance to self!
  v1.mat <- replicate(length(v2.dist), v1.dist)
  v2.mat <- t(replicate(length(v1.dist), v2.dist))
  return((v1.mat + v2.mat) + path.length)
}
