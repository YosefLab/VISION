#' Cluster single cells so signal is maintained but the sample size and noise
#' are reduce, using the louvain clustering algorithm
#' @param exprData the expression data matrix
#' @param cellsPerPartition control over the minimum number of cells to put into each supercell
#' @param preserve_clusters named factor vector denoted cluster boundaries to preserve
#' @param latentSpace numeric matrix cells x components
#' @importFrom Matrix tcrossprod
#' @importFrom Matrix rowMeans
#' @return pooled cells - named list of vectors - cells in each supercell
applyMicroClustering <- function(
                         exprData, cellsPerPartition=100,
                         filterInput = "fano",
                         filterThreshold = round(ncol(exprData)*0.2),
                         preserve_clusters = NULL,
			 random=F,
                         latentSpace = matrix(NA, 1, 1)) {

    if (all(dim(latentSpace) == c(1, 1))) {
        exprData <- matLog2(exprData)
        gene_passes <- applyFilters(exprData, filterThreshold, filterInput)
        fexpr <- exprData[gene_passes, ]

        # Compute wcov using matrix operations to avoid
        # creating a large dense matrix

        N <- ncol(fexpr)
        wcov <- tcrossprod(fexpr) / N

        mu <- as.matrix(rowMeans(fexpr), ncol = 1)
        mumu <- tcrossprod(mu)
        wcov <- as.matrix(wcov - mumu)

        # SVD of wieghted correlation matrix
        ncomp <- min(ncol(fexpr), nrow(fexpr), 10)
        decomp <- rsvd::rsvd(wcov, k = ncomp)
        evec <- t(decomp$u)

        # Project down using computed eigenvectors
        res <- (evec %*% fexpr) - as.vector(evec %*% mu)
        res <- as.matrix(res)
        res <- t(res) # avoid transposing many times below
    } else {
        res <- latentSpace
    }

    # If 'preserve_clusters' is provided, use these as the pre-clustering
    #   Otherwise, compute clusters using knn graph and louvain
    if (!is.null(preserve_clusters)) {
        cluster_names <- levels(preserve_clusters)
        cl <- lapply(cluster_names, function(level){
            cells_at_level <- names(preserve_clusters)[preserve_clusters == level]
            return(cells_at_level)
        })
        names(cl) <- cluster_names

    } else {
        n_workers <- getWorkerCount()
        kn <- ball_tree_knn(res,
                            min(round(sqrt(nrow(res))), 30),
                            n_workers)

        cl <- louvainCluster(kn, res)
    }

    pools <- readjust_clusters(cl, res, cellsPerPartition = cellsPerPartition)

    # Rename clusters
    cn <- paste0("microcluster ", 1:length(pools))
    names(pools) <- cn

    return(pools)

}


#' Aggregate meta-data for cells in pools
#'
#' @param metaData data.frame of meta-data for cells
#' @param pools named list of character vector describing cell ids
#' in each pool
#' @return a data.frame of new meta-data
createPooledMetaData <- function(metaData, pools) {

    poolMetaData <- data.frame(row.names = names(pools))

    for (sigName in colnames(metaData)) {

        scores <- metaData[[sigName]]
        names(scores) <- rownames(metaData)

        if (is.factor(scores)){
            ## Need to compute majority level in each group
            N_MicroClusters <- nrow(poolMetaData)
            newlevels <- union("~", levels(scores))

            # vector to store the final levels in
            clustScores <- factor(integer(N_MicroClusters),
                            levels = newlevels)
            names(clustScores) <- rownames(poolMetaData)

            # vectors for the proportion of each cluster level
            clustScoresLevels <- list()
            for (level in levels(scores)) {
                clustScoresL <- numeric(N_MicroClusters)
                names(clustScoresL) <- rownames(poolMetaData)
                clustScoresLevels[[level]] <- clustScoresL
            }

            for (clust in names(pools)){

                pool <- pools[[clust]]

                if (length(pool) == 0){ #TODO: This shouldn't happen
                    clustScores[clust] <- "~"
                    next
                }

                vals <- scores[match(pool, names(scores))]
                freq <- table(vals) / length(vals)
                maxval <- freq[which.max(freq)]

                if (maxval >= .5){
                    clust_val <- names(maxval)
                } else {
                    clust_val <- "~"
                }
                clustScores[clust] <- clust_val

                for (level in names(freq)){
                    clustScoresLevels[[level]][clust] <- freq[level]
                }
            }

            poolMetaData[[sigName]] <- clustScores

            for (level in names(clustScoresLevels)){

                newName <- paste(sigName, level, sep = "_")
                poolMetaData[[newName]] <- clustScoresLevels[[level]]

            }
        } else { # Then it must be numeric, just average

            ## Need to compute majority level in each group
            N_MicroClusters <- nrow(poolMetaData)

            clustScores <- numeric(N_MicroClusters)
            names(clustScores) <- rownames(poolMetaData)

            for (clust in names(clustScores)){

                pool <- pools[[clust]]

                if (length(pool) == 0){ #TODO: This shouldn't happen
                    clustScores[clust] <- 0
                    next
                }

                vals <- scores[match(pool, names(scores))]
                clust_val <- mean(vals)
                clustScores[clust] <- clust_val
            }

            poolMetaData[[sigName]] <- clustScores
        }
    }

    return(poolMetaData)

}

#' Applies the Louvain algorithm to generate micro-clustered data
#'
#' @param kn List of nearest neighbor indices and euclidean distances to these
#' nearest neighbors
#' @param data Data matrix
#' @return List of clusters, each entry being a vector of indices representing
#' samples in the cluster.
louvainCluster <- function(kn, data) {

    nn <- kn[[1]]
    d <- kn[[2]]
    sigma <- apply(d, 1, max)
    d <- exp(-1 * (d*d) / sigma^2)

    nnl <- lapply(1:nrow(nn), function(i) nn[i,])

    # Create an undirected knn graph
    g <- igraph::graph_from_adj_list(nnl, mode="out")
    igraph::E(g)$weights <- as.vector(t(d))
    g <- igraph::as.undirected(g, mode="each")

    # Now apply the louvain algorithm to cluster the graph
    cl <- igraph::cluster_louvain(g)

    # Gather cluster vector to list of clusters
    clusters <- list()
    mem <- as.vector(igraph::membership(cl))
    for (i in 1:length(mem)) {
        n <- as.character(mem[[i]])
        if (n %in% names(clusters)) {
            clusters[[n]] <- c(clusters[[n]], i)
        } else {
            clusters[[n]] <- c(i)
        }
    }

    clusters <- lapply(clusters, function(i) i <- rownames(data)[i])

    return(clusters)

}

#' Repartitions existing clusters to achieve desired granularity.
#'
#' By default, minimum number of clusters to be generated is the squareroot of
#' the number of cells.
#' @importFrom stats kmeans
#' @param clusters List of clusters, each entry being a vector of cells in a
#' cluster.
#' @param data NUM_SAMPLES x NUM_FEATURES data matrix that was used to generate
#' clusters
#' @param cellsPerPartition the number of cells for a single partition the
#' algorithm should aim for
#' @return Repartitioned clusters, such that a desireable number of
#' microclusters is acheived.
readjust_clusters <- function(clusters, data, cellsPerPartition=100) {

    NUM_PARTITIONS = round(nrow(data) / cellsPerPartition)
    EPSILON = .15

    currPart = length(clusters)
    clusterList <- list()

    while (currPart < ((1 - EPSILON)*NUM_PARTITIONS)) {
        clusterList <- list()
        cluster_offset = 0
        for (i in 1:length(clusters)) {

            # Apply kmeans clustering to existing cluster
            currCl = clusters[[i]]
            subData <- data[currCl,]
            if (length(currCl) > cellsPerPartition) {
                nCl <- kmeans(subData,
                              centers=round(nrow(subData) / cellsPerPartition),
                              iter.max=100)
            } else {
                nCl <- kmeans(subData, centers=1, iter.max=100)
            }
            newClust <- nCl$cluster
            # Gather cluster vector to list of clusters
            for (i in 1:length(newClust)) {
                n <- as.character(newClust[[i]] + cluster_offset)
                sample_n <- names(newClust)[[i]]
                if (n %in% names(clusterList)) {
                    clusterList[[n]] <- c(clusterList[[n]], sample_n)
                } else {
                    clusterList[[n]] <- c(sample_n)
                }

            }

            # Now add to cluster offset for next re-clustering
            cluster_offset <- cluster_offset + max(newClust)
        }

        currPart <- length(clusterList)
        clusters <- clusterList

    }

    return(clusters)
}


#' create "super-cells" by pooling together single cells
#'
#' This function calls createPools on batches of clusters in parallel
#'
#' @param cl cluster association of each cell
#' @param expr expression data
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @return a matrx of expression data for the pooled cells
createPoolsBatch <- function(cl, expr) {

    availableCores <- max(parallel::detectCores() - 1, 1)
    n_workers <- min(availableCores, 10)

    cl_batches <- batchify(cl, 500, n_workers = n_workers)

    pool_batches <- parallel::mclapply(cl_batches, function(cl_batch) {
                       pool_expression <- createPools(cl_batch, expr)
                       return(pool_expression)
    }, mc.cores = min(n_workers, length(cl_batches)))

    pool_data <- do.call(cbind, pool_batches)

    return(pool_data)
}

#' create "super-cells" by pooling together single cells
#' @param cl cluster association of each cell (list of vector of column names)
#' @param expr expression data (genes x cells matrix)
#' @return a matrx of expression data for the pooled cells (genes x pools)
createPools <- function(cl, expr) {

    # Need to construct a large cells x pools matrix
    cl_data <- lapply(seq(length(cl)), function(i) {
        cluster <- cl[[i]]
        cell_indices <- match(cluster, colnames(expr))
        pool_indices <- rep(i, length(cell_indices))
        return(list(cell_indices, pool_indices))
    })

    i <- unlist(lapply(cl_data, function(x) x[[1]]))
    j <- unlist(lapply(cl_data, function(x) x[[2]]))

    dimnames <- list(
                     colnames(expr),
                     names(cl)
                     )
    dims <- c(ncol(expr), length(cl))

    poolSparseMatrix <- sparseMatrix(i = i, j = j,
                                    dims = dims,
                                    dimnames = dimnames)
    pool_data <- as.matrix(expr %*% poolSparseMatrix)
    cl_sizes <- vapply(cl, length, FUN.VALUE = 1)
    pool_data <- t(t(pool_data) / cl_sizes)  # BAH, maybe use sparse values here?

    return(pool_data)
}
