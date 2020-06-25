#' Pool cells into microclusters
#'
#' Merges similar transcriptional profiles into representative 'pools'
#'
#' A latent space is computed for the expression data via PCA after
#' filtering on genes (using parameters \code{filterInput} and \code{filterThreshold}).
#'
#' Alternately, a latent space can be supplied via the \code{latentSpace} argument
#'
#' Euclidean distance within the latent space is then used to create cell pools
#'
#'
#' @param exprData the expression data matrix
#' @param cellsPerPartition control over the target number of cells to put into each supercell
#' @param filterInput name of filtering method ('threshold' or 'fano') or list of
#' genes to use when computing projections.
#' @param filterThreshold Threshold to apply when using the 'threshold' or 'fano' projection genes filter.
#' If greater than 1, this specifies the number of cells in which a gene must be detected
#' for it to be used when computing PCA. If less than 1, this instead specifies the proportion of cells needed
#' @param filterNumMad Number of median absolute deviations to use when selecting highly-variable
#' genes in each mean-sorted bin of genes
#' @param latentSpace (Optional) Latent space to be used instead of PCA numeric matrix cells x components
#' @param K Number of neighbors to use for finding pools.
#' @importFrom Matrix tcrossprod
#' @importFrom Matrix rowMeans
#' @return pooled cells - named list of vectors - cells in each supercell
#' @export
applyMicroClustering <- function(
                         exprData, cellsPerPartition=10,
                         filterInput = "fano",
                         filterThreshold = round(ncol(exprData) * 0.05),
                         filterNumMad = 2,
                         latentSpace = NULL, K=round(sqrt(ncol(exprData)))) {

    if (is.data.frame(exprData)){
        exprData <- data.matrix(exprData)
    }

    if (is.null(latentSpace) || all(dim(latentSpace) == c(1, 1))) {
        exprData <- matLog2(exprData)


        message("    Computing a latent space for microclustering using PCA...")
        if (length(filterInput) > 1){
            gene_passes <- intersect(filterInput, rownames(exprData))
            if (length(gene_passes) == 0){
                stop("Supplied list of genes in `filterInput` does not match any rows of `exprData`")
            } else {
                message(
                    sprintf("    Using supplied list of genes: Found %i/%i matches", length(gene_passes), length(filterInput))
                    )
            }
        } else {
            message("    Determining lateng space genes...")
            gene_passes <- applyFilters(exprData, filterInput, filterThreshold, filterNumMad)

            if (length(gene_passes) == 0){
                stop(
                    sprintf("Filtering with (filterInput=\"%s\", filterThreshold=%i) results in 0 genes\n  Set a lower threshold and re-run", filterInput, filterThreshold)
                    )
            }
        }

        fexpr <- exprData[gene_passes, , drop = FALSE]

        # Compute wcov using matrix operations to avoid
        # creating a large dense matrix

        message("    Performing PCA...")
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
        message(
            sprintf("    Using supplied latent space with %i components", ncol(res))
            )
    }

    message("    Performing initial coarse-clustering...")

    kn <- find_knn_parallel(res, min(K, 30))

    cl <- louvainCluster(kn, res)

    message("    Further partitioning coarse clusters...")
    pools <- readjust_clusters(cl, res, cellsPerPartition = cellsPerPartition)

    # Rename clusters
    cn <- paste0("microcluster_", 1:length(pools))
    names(pools) <- cn

    message(
        sprintf("    Micro-pooling completed reducing %i cells into %i pools",
            nrow(res), length(pools))
        )
    return(pools)

}


#' Aggregate meta-data for cells in pools
#'
#' Used to pool a meta-data data.frame which may contain a mixture
#' of numeric and factor variables
#'
#' For numerical variables, the pooled value is just the average of
#' the value for cells in the pool
#'
#' For factors, the pooled factor value represents the majority level
#' across cells in the pool.  If there is no simple majority (e.g. no
#' factor level > 50\%), then the level '~' is substituted.
#'
#' Additionally, for factor variables, a new variable of the form
#' \code{<Variable>_<Level>} is created for each level in the factor with a value
#' equal to the proportion of cells in the pool with that factor level.
#'
#' @param metaData data.frame of meta-data for cells
#' @param pools named list of character vector describing cell ids
#' in each pool
#' @return A data.frame of pooled meta-data
#' @export
poolMetaData <- function(metaData, pools) {

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
    sigma <- apply(d, 1, function(x) quantile(x, c(.5))[[1]])
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
#' @param data NUM_SAMPLES x NUM_PROTEINS data matrix that was used to generate
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

#' Pools columns of a numeric matrix
#'
#' Uses the provided pools to merge columns of the supplied data matrix
#'
#' This would typically be used on a gene expression matrix (genes X cells) to
#' pool cells.
#'
#' The \code{pools} argument is obtained by running \code{applyMicroClustering}
#'
#' @param data data.frame or matrix to pool
#' @param pools named list of character vector describing cell ids
#' in each pool
#' @importFrom parallel mclapply
#' @return Matrix of pooled data
#' @export
poolMatrixCols <- function(data, pools) {

    n_workers <- getOption("mc.cores")
    n_workers <- if (is.null(n_workers)) 2 else n_workers

    cl_batches <- batchify(pools, 500, n_workers = n_workers)

    pool_batches <- parallel::mclapply(cl_batches, function(cl_batch) {
                       pooled_batch <- poolMatrixCols_Inner(data, cl_batch)
                       return(pooled_batch)
    })

    pooled_data <- do.call(cbind, pool_batches)

    return(pooled_data)
}


#' Pools rows of a numeric matrix
#'
#' Uses the provided pools to merge rows of the supplied data matrix
#'
#' Same as poolMatrixCols only operates on rows.
#'
#' The \code{pools} argument is obtained by running \code{applyMicroClustering}
#'
#' @param data data.frame or matrix to pool
#' @param pools named list of character vector describing cell ids
#' in each pool
#' @return Matrix of pooled data
#' @export
poolMatrixRows <- function(data, pools) {

    data <- t(data)

    pooled_data <- poolMatrixCols(data, pools)

    pooled_data <- t(pooled_data)

    return(pooled_data)
}

#' create "super-cells" by pooling together single cells
#' @param expr expression data (genes x cells matrix)
#' @param pools cluster association of each cell
#' @return a matrx of expression data for the pooled cells (genes x pools)
poolMatrixCols_Inner <- function(expr, pools) {

    # Need to construct a large cells x pools matrix
    cl_data <- lapply(seq(length(pools)), function(i) {
        cluster <- pools[[i]]
        cell_indices <- match(cluster, colnames(expr))
        pool_indices <- rep(i, length(cell_indices))
        return(list(cell_indices, pool_indices))
    })

    i <- unlist(lapply(cl_data, function(x) x[[1]]))
    j <- unlist(lapply(cl_data, function(x) x[[2]]))

    dimnames <- list(
                     colnames(expr),
                     names(pools)
                     )
    dims <- c(ncol(expr), length(pools))

    poolSparseMatrix <- sparseMatrix(i = i, j = j,
                                    dims = dims,
                                    dimnames = dimnames)
    pool_data <- as.matrix(expr %*% poolSparseMatrix)
    cl_sizes <- vapply(pools, length, FUN.VALUE = 1)
    pool_data <- t(t(pool_data) / cl_sizes)

    return(pool_data)
}


#' Performs a binary search on a depth d such that 
#' if depth(u, v) <= d then u and v are in the same cluster
#'
#' @param tree object of class phylo
#' @param reach number of clusters to attempt to generate
#' @return List of clusters, each entry being a vector of indices representing
#' samples in the cluster.
treeCluster <- function(tree, reach=10) {
    high <- length(tree$tip.label)
    low <- 0
    while (T) {
        if (high == low) {
            break
        }
        d <- round((high + low) / 2)
        if (d == 0) {
          break
        }
        
        num_clusters <- 0
        seen_ancestors <- c()
        cl <- list()

        for (cell in seq_len(length(tree$tip.label))) {
            
            ancestor <- ancestor_at_depth(tree, cell, d)
            cell <- tree$tip.label[cell]
            if (ancestor %in% seen_ancestors) {
              cl[[match(ancestor, seen_ancestors)]] <- append(cl[[match(ancestor, seen_ancestors)]], cell)
            } else {
                num_clusters <- num_clusters + 1
                seen_ancestors <- append(seen_ancestors, ancestor)
                cl[[num_clusters]] <- c(cell)
            }
        }

        if (num_clusters >= reach) {
            if(low == d) {
                break
            }
            low <- d
        } else if (num_clusters < reach) {
            if(high == d) {
                break
            }
            high <- d
        }
    }

    return(cl)
}


#' Performs a binary search on a depth d such that 
#' if depth(u, v) <= d then u and v are in the same cluster
#'
#' @param tree object of class phylo
#' @param reach number of clusters to attempt to generate
#' @return List of clusters, each entry being a vector of indices representing
#' samples in the cluster.
treeCluster2 <- function(tree, reach=10) {
    if (reach > length(tree$tip.label)) {
        stop("Number of clusters is too high.")
    }
    
    node_depths <- node.depth(tree)
    root <- find_root(tree)
    cluster_parents <- c()
    cluster_parents[[as.name(root)]] <- node_depths[root]
    
    # get the top level internal nodes
    while (T) {
      cluster_parents <- cluster_parents[order(unlist(cluster_parents), decreasing = T)]
      remove <- as.integer(names(cluster_parents)[1])
      cluster_parents <- cluster_parents[-1]
      
      children <- get_children(tree, remove)
      for (child in children) {
          cluster_parents[[as.name(child)]] <- node_depths[child]
      }
      
      if (length(cluster_parents) >= reach) {
          break
      }
    }
    
    cl <- list()
    for (cluster in seq_len(length(cluster_parents))) {
        cellId <- as.integer(names(cluster_parents)[cluster])
      
        all_children <- get_all_children(tree, cellId) %>% (function(x) {return(tree$tip.label[x])})
        cl[[cluster]] <- all_children
    }
    
    return(cl)
}