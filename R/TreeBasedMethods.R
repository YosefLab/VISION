#' Applies the Simple PPT algorithm onto the expression data.
#'
#' @param exprData Expression data -- Num_Genes x Num_Samples
#' @param numCores Number of cores to use during this analysis
#' @param permExprData a list of permutated expression datasets,
#' to use for significance estimation of the tree [default:NULL]
#' @param nNodes_ Number of nodes to find. Default is sqrt(N)
#' @param sigma regularization parameter for soft-assignment of data points to
#' nodes, used as the variance
#'          of a guassian kernel. If 0, this is estimated automatically
#' @param gamma graph-level regularization parameter, controlling the tradeoff
#' between the noise-levels
#'          in the data and the graph smoothness. If 0, this is estimated
#'          automatically.
#' @return Information on the fitten tree
#'      \itemize{
#'     \item C: spatial positions of the tree nodes in NUM_FEATURES dimensional
#'     space
#'     \item W: Unweighted (binary) adjacency matrix of the fitten tree
#'     \item distMat: distance matrix between each tree node to each datapoint
#'     \item mse: the Mean-Squared-Error of the fitten tree
#'     \item zscore: a significance score for the fitted tree
#'      }

applySimplePPT <- function(exprData, numCores, permExprData = NULL,
                            nNodes_ = round(sqrt(ncol(exprData))), sigma=0, gamma=0) {
    MIN_GAMMA <- 1e-5
    MAX_GAMMA <- 1e5
    DEF_TOL <- 1e-2
    DEF_MAX_ITER <- 50

    C <- NULL
    Wt <- NULL

    if (sigma == 0) {
    km <- stats::kmeans(t(exprData), centers=round(sqrt(ncol(exprData))),
                        nstart=1, iter.max=50)$centers

    sigma <- mean(apply(as.matrix(sqdist(t(exprData), km)), 1, min))
    }

    if (gamma == 0) {

    currGamma <- MIN_GAMMA
    nNodes <- round(log(ncol(exprData)))

    prevMSE <- -Inf
    minMSE <- Inf
    minMSEGamma <- MIN_GAMMA

    tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
    C <- tr$C
    Wt <- tr$W
    currMSE <- tr$mse
    minGamma <- MIN_GAMMA

    while ( ((prevMSE / currMSE) - 1 < 0.05) && currGamma <= MAX_GAMMA) {
        prevMSE <- currMSE
        currGamma <- currGamma * 10
        tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
        C <- tr$C
        Wt <- tr$W
        currMSE <- tr$mse
        if (currMSE < minMSE) {
        minMSE <- currMSE
        minMSEGamma <- currGamma
        }
    }

    if (currGamma == MAX_GAMMA) {
        currGamma <- minGamma
        tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
    }
    minGamma <- currGamma

    while( (currMSE < prevMSE) && ((prevMSE / currMSE) - 1 > 0.05) ) {
        prevMSE <- currMSE
        currGamma <- currGamma * 10
        minGamma <- minGamma * (10^(1/3))
        tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
        C <- tr$C
        Wt <- tr$W
        currMSE <- tr$mse
    }

    if (nNodes_ > nNodes) {

        if (nNodes_ != 0) {
        nNodes <- nNodes_
        } else {
        nNodes <- round(sqrt(ncol(exprData)))
        }

        tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
        C <- tr$C
        Wt <- tr$W
        currMSE <- tr$mse

        # Calculate the degree distribution
        deg <- colSums(Wt)
        degDist <- unname(table(deg))
        deg_g2c <- 0
        if (length(degDist) > 2) {
        deg_g2c <- sum(degDist[seq(3, length(degDist))])
        }
        deg_g2f <- deg_g2c / nNodes

        while ( !(deg_g2c > 0 && (deg_g2f <= 0.1 || deg_g2c < 5)) && (currGamma >= minGamma) ) {
        currGamma <- currGamma / sqrt(10)
        tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
        C <- tr$C
        Wt <- tr$W
        currMSE <- tr$mse

        deg <- colSums(Wt)
        degDist <- unname(table(deg))
        deg_g2c <- 0
        if (length(degDist) > 2) {
            deg_g2c <- sum(degDist[seq(3, length(degDist))])
        }
        deg_g2f <- deg_g2c / nNodes

        }

    }

    gamma <- currGamma

    }

    tr <- fitTree(exprData, nNodes, sigma, gamma, DEF_TOL, DEF_MAX_ITER)

    if (!is.null(permExprData)) {
    mses <- sapply(permExprData, function(permdata) {
        return(fitTree(permdata, nNodes, sigma, gamma, DEF_TOL, DEF_MAX_ITER)$mse)
    })
    zscore <- log1p((tr$mse - mean(mses)) / stats::sd(mses))
    } else {
    zscore <- NULL
    }


    return(list(princPnts = tr$C, adjMat = tr$W, distMat = sqdist(t(exprData), t(C)),
                mse = tr$mse, zscore = zscore))
}

#' Fit tree using input parameters
#'
#' @importFrom matrixStats logSumExp
#' @param expr Data to fit (NUM_GENES x NUM_SAMPLES)
#' @param nNodes Number of nodes in the fitted tree, default is square-root of number of data points
#' @param sigma Regularization parameter for soft-assignment of data points to nodes, used as the
#'                  variance of a gaussian kernel. If 0, this is estimated automatically.
#' @param gamma Graph-level regularization parameter, controlling the tradeoff between the noise-levels
#'                  in the data and the graph smoothness. If 0, this is estimated automatically
#' @param tol Tolerance to use when fitting the tree
#' @param maxIter Maximum number of Iterations ot run the algorithm for
#' @return (list) Tree characteristics:
#'      \itemize{
#'     \item C spatial positions of the tree nodes in NUM_FEATURES dimensional space
#'     \item unweighted (binary) adjacency matrix
#'     \item the mean-squared error of the tree
#'      }

fitTree <- function(expr, nNodes, sigma, gamma, tol, maxIter) {
    km <- stats::kmeans(t(expr), centers=nNodes, nstart=10, iter.max=100)$centers
    cc_dist <- as.matrix(sqdist(km, km))
    cx_dist <- as.matrix(sqdist(t(expr), km))
    prevScore = Inf
    currScore = -Inf
    currIter = 0

    while (!(prevScore - currScore < tol) && !(currIter > maxIter)) {
    currIter <- currIter + 1
    prevScore <- currScore
    W <- igraph::mst(igraph::graph_from_adjacency_matrix(cc_dist,
                                                            weighted= TRUE,
                                                            mode="undirected"))
    Wt <- igraph::get.adjacency(W, sparse=FALSE)

    Ptmp <- -(cx_dist / sigma)
    Psums <- matrix(rep(apply(Ptmp, 1, logSumExp), each=ncol(Ptmp)),
                    nrow=nrow(Ptmp), ncol=ncol(Ptmp), byrow=TRUE)
    P <- exp(Ptmp - Psums)

    delta <- diag(colSums(P))
    L <- igraph::laplacian_matrix(W)
    xp <- crossprod(t(expr), P)
    invg <- as.matrix(solve( ((2 / gamma) * L) + delta))
    C <- tcrossprod(xp, invg)

    cc_dist <- as.matrix(dist.matrix(t(C)))
    cx_dist <- as.matrix(sqdist(t(expr), t(C)))

    P <- clipBottom(P, mi=min(P[P>0]))
    currScore <- sum(Wt * cc_dist) + (gamma * sum(P * ((cx_dist) + (sigma * log(P)))))

    }

    return(list(C = C, W = Wt, mse = getMSE(C, expr)))

}



#' Calculates the MSE between C and X
#'
#' @param C d x m matrix
#' @param X d x n matrix
#' @return Mean Squared Error between C and X.
getMSE <- function(C, X) {

    if (is.na(C) || is.na(X)) {
    return(NULL)
    }
    mse <- mean( apply( as.matrix(sqdist(t(X), t(C))), 1, min))
    return(mse)
}



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
#' @return (list) projection information:
#'      \itemize{
#'     \item{"spatial"}{The D-dimensional position of the projected data points}
#'     \item{"edge"}{a Nx2 matrix, were line i has the indices identifying the edge that datapoint i was projected on,
#'     represented as (node a, node b). For consistency and convenience, it is maintained that a < b}
#'     \item{"edgePos"}{an N-length numeric with values in [0,1], the relative position on the edge of the datapoint.
#'     0 is node a, 1 is node b, .5 is the exact middle of the edge, etc.}
#'      }
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
                edgePos = projections[3,],
                spatial = projections[-c(1:3),]))
}



#' Calculate distance matrix between all pairs of ponts based on their projection onto the tree
#'
#' @param princPnts (D x K numeric) the spatial locations of the principle points (tree vertices)
#' @param princAdj (K x K logical) adjacency matrix of the principle graph
#' @param edgeAssoc (2 x N) for each point, the edge it is projected to (represented as (V1,V2), where V1<V2)
#' @param edgePos (length N, numeric) relative postion on the edge for each point, in range [0,1]
#'
#' @return non-negative symmetric matrix in which [i,j] is the tree-based distance between points i, j.
calculateTreeDistances <- function(princPnts, princAdj, edgeAssoc, edgePos) {
    # get all distances in principle tree
    princAdjW <- sqdist(t(princPnts), t(princPnts)) * princAdj

    princGraph <- igraph::graph_from_adjacency_matrix(princAdjW,
                                                    weighted = TRUE,
                                                    mode = "undirected")
    nodeDistmat <- igraph::distances(princGraph)

    princEdges <- apply(igraph::get.edgelist(princGraph), 1, as.numeric)
    edgeToPnts <- apply(princEdges, 2, function(x) { apply(edgeAssoc==x, 2, all) })

    distmat <- matrix(rep(NA, NROW(edgeToPnts) ^ 2), NROW(edgeToPnts))
    ## compute intra-edge distances. Store in list for later calclations and set
    ## values in result matrix
    ## loop contains assignment to external matrix, so apply can't be used.
    ## Alternative is to use Reduce on lapply result, but memory footprint could
    ## be problematic (it's K-choose-2 matrices of size NxN)
    intraDist <- list()
    for (i in 1:NCOL(princEdges)) {
    inEdgeDist <- calcIntraEdgeDistMat(edge.len = nodeDistmat[princEdges[1,i],
                                                                princEdges[2,i]],
                                        edgePos = edgePos[edgeToPnts[,i]])
    intraDist[[i]] <- inEdgeDist[,c(1,NCOL(inEdgeDist))]
    distmat[edgeToPnts[,i], edgeToPnts[,i]] <- inEdgeDist[-c(1,NROW(inEdgeDist)),
                                                            -c(1,NROW(inEdgeDist))]
    }

    ## for each pair of edges, calculate inter-edge distances and set them in the empty matrix
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
    return(pmax(distmat, t(distmat), na.rm = TRUE))
}

#' Calculate distances between all points on a given edge, including edge vertices
#' @param edge.len the length of the node
#' @param edgePos the relative positions of the points on the edge (in range [0,1]).
#' @return a distance matrix, where the first and last points are the edge vertices
calcIntraEdgeDistMat <- function(edge.len, edgePos) {
    edgePos <- c(0, edgePos, 1) * edge.len
    pos.mat <- replicate(length(edgePos), edgePos)
    return(abs(pos.mat - t(pos.mat)))
}

#' Calculate all distances between points on two different edges
#' @param v1.dist a vectors of all distances between all points on the first
#' edge, and the vertex of the first edge that is closer to the second edge
#' @param v2.dist a vectors of all distances between all points on the second
#' edge, and the vertex of the second edge that is closer to the first edge
#' @param path.length the length of the path connected the two vertices
#' @return a distamce matrix between all points on edge1 and all points on edge2
calcInterEdgeDistMat <- function(v1.dist, v2.dist, path.length) {
    # note that input vector must not contain distance to self!
    v1.mat <- v1.dist %o% rep(1, length(v2.dist))
    v2.mat <- t(v2.dist %o% rep(1, length(v1.dist)))
    return((v1.mat + v2.mat) + path.length)
}

#' Find K nearest neighbors for each vector in query among the vectors in the
#' given data
#' @param data the search-space to look for nearest neighbors
#' @param query the points to find neighbors for (not necessarily in the data)
#' @param k the number of neighbors to find for each vectors
#' @param BPPARAM the parallelization backend to use
#' @return a list of the k neughbors for each of the vectors in `query`
findNeighbors <- function(data, query, k, BPPARAM=bpparam()) {

    neighborhood <- lapply(1:ncol(query), function(x) {
    vkn <- ball_tree_vector_knn(t(data), query[,x], k, BPPARAM$workers)
    return(vkn)
    })
    return(neighborhood)
}
