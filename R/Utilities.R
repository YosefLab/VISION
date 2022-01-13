#' Utility methods


#' Helper utility to group list items into batches
#'
#' This is used inside the batchSigEval function
#'
#' @param items list of items to group into batches
#' @param per_batch approximate target for size of batches
#' @param n_workers number of batches is made divisible by n_workers
#' @return list of list of items
batchify <- function(items, per_batch, n_workers = 1) {
    if (length(items) == 0){
        return(list(list()))
    }
    n_iterations <- round(length(items) / (per_batch * n_workers))

    if (n_iterations == 0){
        n_batches <- ceiling(length(items) / per_batch)
    } else {
        n_batches <- n_iterations * n_workers
    }
    per_batch <- ceiling(length(items) / n_batches)
    n_batches <- ceiling(length(items) / per_batch)

    out <- lapply(seq_len(n_batches), function(i) {
        start_i <- (i - 1) * per_batch + 1
        end_i <- i * per_batch
        if (end_i > length(items)){
            end_i  <- length(items)
        }
        return(items[start_i:end_i])
    })

    return(out)
}


#' Checks the version of the Vision object and displays error if necessary
#'
#' @param object Vision object
#' @return NULL
versionCheck <- function(object) {

    templateStr <- paste0(
        "This Vision object was created with an older version of the library.",
        "  To view, either install an older version (commit #COMMIT) from GitHub (www.github.com/YosefLab/Vision) or",
        " recreate the object and re-run analyze"
    )

    if (!.hasSlot(object, "version")) {
        msg <- gsub("#COMMIT", "0cd5268", templateStr)
        stop(msg, call. = FALSE)
    }

    if(object@version < 1.0) {
        msg <- gsub("#COMMIT", "0cd5268", templateStr)
        stop(msg, call. = FALSE)
    }

    if(object@version < 1.1) {
        msg <- gsub("#COMMIT", "2db4552", templateStr)
        stop(msg, call. = FALSE)
    }

    if(object@version < 1.2) {
        msg <- gsub("#COMMIT", "5c085cb", templateStr)
        stop(msg, call. = FALSE)
    }

    # Add new commit hashes here as version increases and breaks backward compatibility

    return()
}


#' log2-scale transform a dense OR sparse matrix
#'
#' This avoids the creation of a dense intermediate matrix when
#' operating on sparse matrices
#'
#' Either performs result <- log2(spmat+1) or if scale = TRUE
#' returns result <- log2(spmat/colSums(spmat)*scaleFactor + 1)
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @param spmat sparse Matrix
#' @param scale boolean - whether or not to scale the columns to sum to `scale_factor`
#' @param scaleFactor if scale = TRUE, columns are scaled to sum to this number
#' @return logmat sparse Matrix
matLog2 <- function(spmat, scale = FALSE, scaleFactor = 1e6) {


    if (scale == TRUE) {
        spmat <- t( t(spmat) / colSums(spmat)) * scaleFactor
    }

    if (is(spmat, "sparseMatrix")) {
        matsum <- summary(spmat)

        logx <- log2(matsum$x + 1)

        logmat <- sparseMatrix(i = matsum$i, j = matsum$j,
                               x = logx, dims = dim(spmat),
                               dimnames = dimnames(spmat))
    } else {
        logmat <- log2(spmat + 1)
    }


    return(logmat)

}


#' inverse log-scale transform a dense OR sparse matrix
#'
#' This avoids the creation of a dense intermediate matrix when
#' operating on sparse matrices
#'
#' Performs result <- exp(spmat) - 1
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @param spmat sparse Matrix
#' @return logmat sparse Matrix
ilog1p <- function(spmat) {


    if (is(spmat, "sparseMatrix")) {
        matsum <- summary(spmat)

        ilogx <- exp(matsum$x) - 1

        ilogmat <- sparseMatrix(i = matsum$i, j = matsum$j,
                                x = ilogx, dims = dim(spmat),
                                dimnames = dimnames(spmat))
    } else {
        ilogmat <- exp(spmat) - 1
    }


    return(ilogmat)

}

#' Vectorized wilcox rank-sums test
#'
#' Given indices in cluster_ii, compute the ranksums test
#' on every column for values in cluster_ii vs the rest.
#'
#' @importFrom matrixStats colCounts
#' @importFrom matrixStats colMaxs
#' @param ranks matrix of ranks, each column representing a separate variable
#' @param cluster_ii numeric vector - indices denoting the group to be compared
#' @param check_na - whether or not to check for NA values (slows computation)
#' @param check_ties - whether or not to check for ties in ranks (slows computation)
#' @return pval - numeric vector, pvalue for each column
#' @return stat - numeric vector, test statistic (AUC) for each column
matrix_wilcox <- function(ranks, cluster_ii,
                          check_na = FALSE, check_ties = FALSE){

    # handle edge case
    if (ncol(ranks) == 0){
        p <- numeric()
        stat <- numeric()
        return(list(pval = p, stat = stat))
    }

    subset <- ranks[cluster_ii, , drop = FALSE]

    not_cluster_ii = setdiff(seq(nrow(ranks)), cluster_ii)
    subset_2 = ranks[not_cluster_ii, , drop=FALSE]

    if (check_na){
        n1 <- nrow(subset) - colCounts(subset, value = NA)
	n2 <- nrow(subset_2) - colCounts(subset_2, value=NA)
    } else {
        n1 <- length(cluster_ii)
        n2 <- nrow(ranks) - n1
    }


    r1 <- colSums(subset, na.rm = T)

    u1 <- r1 - n1 * (n1 + 1) / 2

    u <- u1
    AUC <- u / (n1 * n2)
    AUC <- ifelse(is.infinite(AUC), .5, AUC) # Edge case n2 == 0
    AUC <- ifelse(is.na(AUC), .5, AUC) # Edge case, n1 == 0

    # u to z-score
    m_u <- (n1 * n2) / 2
    if (check_ties){
        sd_u <- vapply(seq_len(ncol(ranks)), function(i){
                    Y <- ranks[, i]
                    n1 <- if (length(n1) == 1) n1 else n1[i]
                    n2 <- if (length(n2) == 1) n2 else n2[i]
                    has_ties <- length(Y) != length(unique(Y))
                    if (!has_ties){
                        sd_ui <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
                        return(sd_ui)
                    }
                    n <- n1 + n2
                    n_ties <- table(Y)
                    sd_ui <- sqrt(
                                  n1 * n2 / 12 *
                                  (n + 1 -
                                   sum(n_ties ** 3 - n_ties) / n / (n - 1)
                                  )
                             )
                    return(sd_ui)
        }, FUN.VALUE = 0.0)


    } else {
        sd_u <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    }
    z <- (u - m_u) / sd_u
    z <- -1 * abs(z) # ensure negative
    z <- z + .5 / sd_u # continuity correction
    z[sd_u == 0] <- 0  # handle case where n1 or n2 = 0
    p <- pnorm(z) * 2

    p[p > 1] <- 1  # Can happen due to continuity correction

    return(list(pval = p, stat = AUC))
}


#' C++ wilcox rank-sums test
#'
#' Given indices in cluster_num, and cluster_denom, compute the ranksums test
#' on every row for values in cluster_num vs cluster_denom.
#'
#' Ties are handled correctly but values are assumed to contain no NAs
#'
#' @importFrom parallel mclapply
#' @param data matrix of values, each row representing a separate variable
#' @param cluster_num numeric vector - indices denoting the group to be compared
#' @param cluster_denom numeric vector - indices denoting the other group to be compared
#' @param jobs integer specifying number of parallel jobs.  Can be used to override
#' the default mc.cores option for when running this within an already-parallel loop
#' @return pval - numeric vector, pvalue for each row
#' @return stat - numeric vector, test statistic (AUC) for each row
matrix_wilcox_cpp <- function(data, cluster_num, cluster_denom,
    jobs = getOption("mc.cores", 1L)) {

    # Subsetting individual rows is bad with a sparse matrix
    # instead we subset chunks at a time

    if ( is(data, "sparseMatrix")){
        batches <- batchify(as.numeric(seq_len(nrow(data))), 100)
        items <- as.numeric(seq_len(nrow(data)))
        res <- mclapply(batches, function(batch){
            dsub <- as.matrix(data[batch, , drop = FALSE])
            res_b <- lapply(seq_len(nrow(dsub)), function(i){
                uz <- wilcox_subset(dsub[i, ], cluster_num, cluster_denom)
                return(unlist(uz))
            })
            return(do.call(rbind, res_b))
        }, mc.cores = jobs)
    } else {
        res <- mclapply(seq_len(nrow(data)), function(i) {
            uz <- wilcox_subset(as.numeric(data[i, ]), cluster_num, cluster_denom)
            return(unlist(uz))
        }, mc.cores = jobs)
    }

    res <- as.data.frame(do.call(rbind, res))
    rownames(res) <- rownames(data)

    res$AUC <- res$U / (length(cluster_num) * length(cluster_denom))

    res$pval <- pnorm(res$Z) * 2
    res$pval[res$pval > 1] <- 1 # Can happen due to continuity correction
    res$pval[is.na(res$pval)] <- 1 # Can happen when no expression in either group

    return(res)
}


#' Perform 1vAll factor analysis given a factor matrix and group definition
#'
#' Given indices in cluster_ii, compute the chisq test
#' on every column for values in cluster_ii vs the rest.
#'
#' @param factorDF dataframe of factors
#' @param cluster_ii numeric vector - indices denoting the group to be compared
#' @return pval - numeric vector, pvalue for each column
#' @return stat - numeric vector, test statistic (Cramer's V) for each column
matrix_chisq <- function(factorDF, cluster_ii) {

    out <- lapply(colnames(factorDF), function(var){

                    if (!is.factor(factorDF[[var]])){
                        stop("Error: matrix_chisq must be called on a factor dataframe")
                    }

                    values <- factorDF[, var, drop = F]
                    values <- droplevels(values)

                    values[, 2] <- 0
                    values[cluster_ii, 2] <- 1
                    values[, 2] <- as.factor(values[, 2])
                    M <- table(values)
                    if (nrow(M) > 1 && ncol(M) > 1){
                        suppressWarnings(
                            out <- chisq.test(M)
                        )
                        pval <- out$p.value
                        n <- sum(M)
                        V <- sqrt(out$statistic / n /
                            min(nrow(M) - 1, ncol(M) - 1)
                        )
                        stat <- V # Cramer's V
                    } else {
                        pval <- 1.0
                        stat <- 0
                    }
                return(list(pval = pval, stat = stat))
    })
    names(out) <- colnames(factorDF)

    pvals <- vapply(out, function(x) x$pval, FUN.VALUE = 0.0)
    stat <- vapply(out, function(x) x$stat, FUN.VALUE = 0.0)

    return(list(pval = pvals, stat = stat))
}


#' Change Gene Identifiers
#'
#' Changes gene identifiers by aggregating expression measures (sum).  This
#' is mainly useful when changing from Ensembl IDs to Gene Symbols
#'
#' This method is made fast by the use of sparse matrix multiplication
#'
#' i.e.: (newIds x oldIds) x (oldIds x cells) = (newIds x cells)
#'
#' @importFrom Matrix sparseMatrix
#' @param exp expression matrix (genes x cells)
#' @param newIds character vector specifying the new identifer that corresponds
#' with each row of the input \code{exp} matrix
#' @return a matrix in which rows with duplicate 'newIds' have been
#' summed together
#' @export
#' @examples
#'
#' exp <- matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3), nrow=3)
#' colnames(exp) <- c("Cell1", "Cell2", "Cell3")
#' print(exp)
#'
#' newIds <- c("GeneA", "GeneA", "GeneB")
#'
#' result <- convertGeneIds(exp, newIds)
#' print(result)
#'
convertGeneIds <- function(exp, newIds){

    if (length(newIds) != nrow(exp)){
        stop("`newIds` must have same length as number of rows in `exp`")
    }

    unique_symbols <- sort(unique(newIds))

    ens_id <- seq(nrow(exp)) # index of the ensemble id
    unique_id <- match(newIds, unique_symbols)

    aggMat <- sparseMatrix(i = unique_id, j = ens_id,
                           dims = c(length(unique_symbols), nrow(exp)),
                           dimnames = list(unique_symbols, NULL))

    exp_sym <- aggMat %*% exp

    return(exp_sym)
}


#' Read 10x Output
#'
#' Loads 10x output counts and converts expression to gene symbols
#'
#' This version takes in three files as inputs:
#' \enumerate{
#'   \item matrix.mtx
#'   \item genes.tsv
#'   \item barcodes.tsv
#' }
#'
#' These files are found in the output of "cellranger count" in a folder
#' that looks like:
#'
#' \code{outs/filtered_gene_bc_matrices/mm10}
#'
#' though with the name of whichever genome you are using instead of 'mm10'
#'
#' @param expression path to matrix.mtx
#' @param genes path to genes.tsv
#' @param barcodes path to barcodes.tsv
#' @param ensToSymbol bool denoting whether or not to perform label conversion
#' @importFrom Matrix readMM
#' @importFrom utils read.table
#' @return sparse count matrix with appropriate row/column names
#' @export
read_10x <- function(expression, genes, barcodes, ensToSymbol = TRUE){

    counts <- readMM(expression)

    gene_data <- read.table(genes, header=FALSE)
    symbols <- gene_data$V2

    barcodes <- readLines(barcodes)
    colnames(counts) <- barcodes
    rownames(counts) <- gene_data$V1

    if (ensToSymbol){
        counts <- convertGeneIds(counts, symbols)
    }

    return(counts)
}


#' Read 10x HDF5 Output
#'
#' Loads 10x output counts and converts expression to gene symbols
#'
#' This version uses the h5 file produced by "cellranger count"
#'
#' This file is typically in a folder that looks like:
#'
#' \code{outs/filtered_gene_bc_matrices_h5.h5}
#'
#' @param h5_file path to h5 file
#' @param ensToSymbol bool denoting whether or not to perform label conversion
#' @importFrom Matrix sparseMatrix
#' @return Return value depends on whether this is a Cellranger v3 or Cellranger v2 output
#'     If it is a v2 output, return is a sparse count matrix with gene expression values
#'     If it is a v3 output, return value is a list with two entries:
#'         Expression: sparse count matrix with gene expression counts (genes x cells)
#'         Antibody: sparse count matrix with antibody capture counts (cells x antibodies)
#' @export
read_10x_h5 <- function(h5_file, ensToSymbol = TRUE){

    if (!requireNamespace("hdf5r", quietly = TRUE)){
      stop("Package \"hdf5r\" needed to load this data object.  Please install it.",
           call. = FALSE)
    }

    h5 <- hdf5r::H5File$new(h5_file)

    is_v3 <- FALSE

    tryCatch({
        genomes <- names(h5)

        if (length(genomes) > 1){
            stop("The supplied h5 file has multiple genomes.  Loading this is not supported by this function")
        }

        genome <- genomes[1]
        if ("features" %in% names(h5[[genome]])) {
            is_v3 <- TRUE
        }
        },
        finally = {
            h5$close_all()
        }
    )

    if (is_v3){
        return(read_10x_h5_v3(h5_file, ensToSymbol))
    } else {
        return(read_10x_h5_v2(h5_file, ensToSymbol))
    }
}


#' Read 10x HDF5 Output - CellRanger 2.0
#'
#' Loads 10x output counts and converts expression to gene symbols
#'
#' This version uses the h5 file produced by "cellranger count"
#'
#' This file is typically in a folder that looks like:
#'
#' \code{outs/filtered_gene_bc_matrices_h5.h5}
#'
#' @param h5_file path to h5 file
#' @param ensToSymbol bool denoting whether or not to perform label conversion
#' @importFrom Matrix sparseMatrix
#' @return sparse count matrix with appropriate row/column names
#' @export
read_10x_h5_v2 <- function(h5_file, ensToSymbol = TRUE){

    if (!requireNamespace("hdf5r", quietly = TRUE)){
      stop("Package \"hdf5r\" needed to load this data object.  Please install it.",
           call. = FALSE)
    }

    h5 <- hdf5r::H5File$new(h5_file)

    tryCatch({
        genomes <- names(h5)

        if (length(genomes) > 1){
            stop("The supplied h5 file has multiple genomes.  Loading this is not supported by this function")
        }

        genome <- genomes[1]

        data <- h5[[paste0(genome, "/data")]][]
        data <- as.numeric(data)

        indices <- h5[[paste0(genome, "/indices")]][]
        indptr <- h5[[paste0(genome, "/indptr")]][]
        dims <- h5[[paste0(genome, "/shape")]][]
        ensIds <- h5[[paste0(genome, "/genes")]][]
        symbols <- h5[[paste0(genome, "/gene_names")]][]
        barcodes <- h5[[paste0(genome, "/barcodes")]][]

        dimnames <- list(ensIds, barcodes)

        counts <- sparseMatrix(i = indices + 1,
            p = indptr, x = data, dims = dims,
            dimnames = dimnames
        )
        if (ensToSymbol){
            counts <- convertGeneIds(counts, symbols)
        }
        },
    finally = {
        h5$close_all()
    })


    return(counts)
}


#' Read 10x HDF5 Output - CellRanger 3.0
#'
#' Loads 10x output counts and converts expression to gene symbols
#'
#' This version uses the h5 file produced by "cellranger count"
#'
#' This file is typically in a folder that looks like:
#'
#'     \code{outs/filtered_feature_bc_matrices_h5.h5}
#'
#' @param h5_file path to h5 file
#' @param ensToSymbol bool denoting whether or not to perform label conversion
#' @importFrom Matrix sparseMatrix
#' @return a list with two items
#'         Expression: sparse count matrix with gene expression counts
#'         Antibody: sparse count matrix with antibody capture counts
#' @export
read_10x_h5_v3 <- function(h5_file, ensToSymbol = TRUE){

    if (!requireNamespace("hdf5r", quietly = TRUE)){
      stop("Package \"hdf5r\" needed to load this data object.  Please install it.",
           call. = FALSE)
    }

    h5 <- hdf5r::H5File$new(h5_file)

    tryCatch({
        genomes <- names(h5)

        if (length(genomes) > 1){
            stop("The supplied h5 file has multiple genomes.  Loading this is not supported by this function")
        }

        genome <- genomes[1]

        data <- h5[[paste0(genome, "/data")]][]
        data <- as.numeric(data)

        indices <- h5[[paste0(genome, "/indices")]][]
        indptr <- h5[[paste0(genome, "/indptr")]][]
        dims <- h5[[paste0(genome, "/shape")]][]
        barcodes <- h5[[paste0(genome, "/barcodes")]][]

        features <- h5[[paste0(genome, "/features")]]

        features_df <- data.frame(
            id = features[["id"]][],
            name = features[["name"]][],
            feature_type = features[["feature_type"]][],
            genome = features[["genome"]][],
            pattern = features[["pattern"]][],
            read = features[["read"]][],
            sequence = features[["sequence"]][]
        )

        ids <- features_df$id
        symbols <- features_df$name

        dimnames <- list(ids, barcodes)

        counts <- sparseMatrix(i = indices + 1,
            p = indptr, x = data, dims = dims,
            dimnames = dimnames
        )

        counts_expression <- counts[
            features_df$feature_type == "Gene Expression",
        ]

        symbols_expression <- symbols[
            features_df$feature_type == "Gene Expression"
        ]

        counts_antibody <- counts[
            features_df$feature_type == "Antibody Capture",
        ]
        counts_antibody <- t(counts_antibody)

        if (ensToSymbol){
            counts_expression <- convertGeneIds(counts_expression, symbols_expression)
        }
        },
    finally = {
        h5$close_all()
    })

    return(
       list(
            Expression = counts_expression,
            Antibody = counts_antibody
            )
       )
}


#' Tests for Unnormalized Data
#'
#' Determines if the VISION object is storing unnormalized data
#'
#' @param object VISION object
#' @return bool whether or not there is unnormalize data
hasUnnormalizedData <- function(object) {
    if (all(dim(object@unnormalizedData) == 1)){
        return(FALSE)
    }

    return(TRUE)
}


#' Gets parameters with defaults
#'
#' Used to get analysis parameters and supply defaults if needed
#'
#' @param object VISION object
#' @param param name of parameter that is desired
#' @param subparam name of second level parameter that is desired
#' @return bool whether or not there is unnormalize data
getParam <- function(object, param, subparam = NULL) {

    useDefault <- tryCatch({
    }, error = function(err){
        return(TRUE)
    })

    useDefault <- !(param %in% names(object@params))
    if ((!useDefault) && (!is.null(subparam))) {
        useDefault <- !(subparam %in% names(object@params[[param]]))
    }

    if (!useDefault){
        if (is.null(subparam)){
            return(object@params[[param]])
        } else {
            return(object@params[[param]][[subparam]])
        }
    } else {
        res <- switch(param,
            "latentSpace" = {
                switch(subparam,
                    "name" = "Latent Space",
                    stop(paste0("Unrecognized subparameter name - ", subparam))
                )
            },
            "name" = "",
            stop(paste0("Unrecognized parameter name - ", param))
        )
        return(res)
    }

}


#' Compute row-wise variance on matrix without densifying
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix rowSums
#' @importFrom matrixStats rowVars
#'
#' @param x expression matrix
#' @return numeric vector row-wise variance
rowVarsSp <- function(x) {
    if (is(x, "matrix")) {
        out <- matrixStats::rowVars(x)
        names(out) <- rownames(x)
    } else {
        rm <- Matrix::rowMeans(x)
        out <- Matrix::rowSums(x ^ 2)
        out <- out - 2 * Matrix::rowSums(x * rm)
        out <- out + ncol(x) * rm ^ 2
        out <- out / (ncol(x) - 1)
    }
    return(out)
}

#' Compute col-wise variance on matrix without densifying
#'
#' @importFrom Matrix colMeans
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom matrixStats colVars
#'
#' @param x expression matrix
#' @return numeric vector col-wise variance
colVarsSp <- function(x) {
    if (is(x, "matrix")) {
        out <- matrixStats::colVars(x)
        names(out) <- colnames(x)
    } else {
        rm <- Matrix::colMeans(x)
        out <- Matrix::colSums(x ^ 2)
        out <- out - 2 * Matrix::rowSums(t(x) * rm)
        out <- out + nrow(x) * rm ^ 2
        out <- out / (nrow(x) - 1)
    }
    return(out)
}

#' Parallel KNN
#'
#' Computes nearest-neighbor indices and distances
#' in parallel.  Query is over all points and results
#' do not include the self-distance.
#'
#' @importFrom RANN nn2
#' @importFrom pbmcapply pbmclapply
#'
#' @param data a Samples x Dimensions numeric matrix
#' @param K number of neighbors to query
#' @return list with two items:
#'     index: Samples x K matrix of neighbor indices
#'     dist: Samples x K matrix of neighbor distances
find_knn_parallel <- function(data, K) {

    workers <- getOption("mc.cores")
    if (is.null(workers)){
        workers <- 1
    }

    per_batch <- round(nrow(data)/workers)

    batches <- batchify(seq_len(nrow(data)), per_batch, n_workers = workers)

    nn_out <- mclapply(batches, function(rows){
        nd <- RANN::nn2(data, query = data[rows, , drop = FALSE], k = K+1)
        return(nd)
    })

    idx <- do.call(rbind, lapply(nn_out, function(nd) nd$nn.idx))
    dists <- do.call(rbind, lapply(nn_out, function(nd) nd$nn.dists))

    idx <- idx[, -1, drop = FALSE]
    dists <- dists[, -1, drop = FALSE]

    return(list(index=idx, dist=dists))
}

#' Helper KNN Function for Trees
#'
#' Computes nearest-neighbor indices and distances
#' in serial.  Query is over all points and results
#' do not include the self-distance.
#'
#' @param leaves leaves to find the nearest neighbors of
#' @param k number of neighbors to query
#' @param distances matrix of Samples x Samples matrix of cophenetic tree distances
#' @return list with two items:
#'     index: Samples x K matrix of neighbor indices
#'     dist: Samples x K matrix of neighbor distances
knn_tree <- function(leaves, k, distances) {
  n <- length(leaves)
  rd <- list(idx = matrix(, length(leaves), k+1), dists = matrix(, length(leaves), k+1))
  rownames(rd$idx) <- leaves
  rownames(rd$dists) <- leaves
  for (leaf in leaves) {
    subsetDistances <- distances[leaf, ]
    random <- runif(length(subsetDistances)) * 0.1
    sorted <- sort(subsetDistances + random)
    nearest <- sorted[1:(k+1)] # don't include myself in my own nearest neighbors
    rd$idx[leaf, ] <- names(nearest) %>% (function(x){match(x, colnames(distances))})
    rd$dists[leaf, ] <- nearest
  }
  return(rd)
}


#' Parallel KNN for Trees
#'
#' Computes nearest-neighbor indices and distances
#' in parallel.  Query is over all points and results
#' do not include the self-distance.
#'
#' @importFrom pbmcapply pbmclapply
#'
#' @param tree an object of class phylo
#' @param K number of neighbors to query
#' @return list with two items:
#'     index: Samples x K matrix of neighbor indices
#'     dist: Samples x K matrix of neighbor distances
find_knn_parallel_tree <- function(tree, K) {
  
  workers <- getOption("mc.cores")
  if (is.null(workers)){
    workers <- 1
  }
  
  n <- length(tree$tip.label)
  per_batch <- round(n/workers)
  
  batches <- batchify(tree$tip.label, per_batch, n_workers = workers)
  
  # set edge lengths to uniform if not specified in tree
  if (is.null(tree$edge.length)) {
    tree$edge.length = rep(1, length(tree$edge))
  }
  # find the distances
  distances <- cophenetic.phylo(tree)
  
  nn_out <- mclapply(batches, function(leaves){
    return(knn_tree(leaves, K, distances))
  })
  
  idx <- do.call(rbind, lapply(nn_out, function(nd) nd$idx))
  dists <- do.call(rbind, lapply(nn_out, function(nd) nd$dists))
  
  idx <- idx[, -1, drop = FALSE]
  dists <- dists[, -1, drop = FALSE]
  
  return(list(index=idx, dist=dists))
}



#'  Generate neighbors and weights for a tree object, based on LCA.
#' 
#'  @param tree object of class phylo
#'  @param minSize the minimum number of neighbors per node
#'  @return the hotspot object
lcaBasedTreeKNN <- function(tree, minSize=20) {
  tips <- tree$tip.label
  nTips <- length(tips)
  neighbors <- data.frame(t(matrix(seq_len(nTips), ncol = nTips, nrow= nTips)))
  rownames(neighbors) <- tips
  
  
  weights <- data.frame(matrix(0, ncol = nTips, nrow= nTips))
  for (tip in seq_len(nTips)) {
    my_neighbors <- minSizeCladeNeighbors(tree, tip, minSize)
    
    weights[tip, my_neighbors] <- 1
  }
  
  neighbors_no_diag <- data.frame(matrix(ncol = nTips -1, nrow= nTips))
  weights_no_diag <- data.frame(matrix(ncol = nTips -1, nrow= nTips))
  
  for (tip in seq_len(nTips)) {
    neighbors_no_diag[tip, ] <- neighbors[tip, -tip]
    weights_no_diag[tip, ] <- weights[tip, -tip]
  }
  
  rownames(neighbors_no_diag) <- tips
  rownames(weights_no_diag) <- tips
  
  colnames(neighbors_no_diag) <- seq_len(nTips-1) - 1
  colnames(weights_no_diag) <- seq_len(nTips-1) - 1
  return(list("neighbors"=as.matrix(neighbors_no_diag), "weights"=as.matrix(weights_no_diag)))
}

#' Generate an ultrametric tree
#'
#' @param tree an object of class phylo
#' @return a tree with edge distances such that it is ultrametric.
ultrametric_tree <- function(tree) {
    tree$edge.length <- rep(1, length(tree$edge[,1]))
    nodeDepths <- node.depth(tree)
    for (edgeI in seq_len(length(tree$edge[,1]))) {
      edge <- tree$edge[edgeI,]
      me <- edge[2]
      parent <- edge[1]
      tree$edge.length[edgeI] <- nodeDepths[parent] - nodeDepths[me]
    }
    return(tree)
}

#' Depth first traversal of a tree starting a specified node
#' 
#' @param tree an rooted tree object of class `phylo`
#' @param postorder return nodes in postorder (i.e., starting from leaves)
#' @param source source from which to start traversal. Defaults to root
depthFirstTraverse <- function(tree, postorder=TRUE, source=NULL) {

    if (is.null(source)) {
        source <- find_root(tree)
    }
    queue <- c(source)
    traversal <- c()
    while (length(queue) > 0) {
        node <- queue[[1]]
        if (!(node %in% traversal)) {
            traversal <- c(node, traversal)
            children <- get_children(tree, node)
            queue <- c(queue, children)
        }
        queue <- queue[-1] # pop off first entry in queue
    }

    if (!postorder) {
        traversal <- rev(traversal)
    }
    return(traversal)

}

#' Find the ancestor of a node above a specific depth
#'
#' @param tree an object of class phylo
#' @param node the index of the node
#' @param depth the depth to search to
#' @return the index of the ancestor node
ancestor_at_depth <- function(tree, node, depth) {
    nodeDepths <- node.depth(tree)
    edges <- tree$edge
    if (depth > max(nodeDepths)) {
        stop("Can't search for a depth greater than max depth of tree")
    }
    while (nodeDepths[node] < depth) {
        node <- edges[, 1][edges[, 2] == node]
    }
    return(node)
}

#' Find the root node
#'
#' @param tree an object of class phylo
#' @return the index of the root node
find_root <- function(tree) {
    return(ancestor_at_depth(tree, 1, length(tree$tip.label)))
}


#' Find the children of a node
#'
#' @param tree an object of class phylo
#' @param node the node to get the children of
#' @return the children
get_children <- function(tree, node) {
    edges <- tree$edge
    children <- edges[, 2][edges[, 1] == node]
    if (length(children) == 0) {
        children <- node
    }
    return(children)
}

#' Get the parent of a node
#' @param tree an object of class phylo
#' @param node the node to get parent
#' 
#' @return the immediate parent of the mode, or the node if it is root
get_parent <- function(tree, node) {
  edges <- tree$edge
  parent <- edges[, 1][edges[, 2] == node]
  if (length(parent) > 0){
    return(parent)
  }
  return(node)
}


#' Get's the nearest >= min size neighbors of a node based on clade structure
#' @param tree an object of class phylo
#' @param tip the tip to find the neighbors of
#' @param minSize the minimum number of neighbors of the node (excludes self)
#' 
#' @return the neighbors
minSizeCladeNeighbors <- function(tree, tip, minSize=20) {
  node <- tip
  while (T) {
    neighbors <- get_all_children(tree, node)
    clade_size <- length(neighbors)
    if (clade_size > minSize) {
      return(neighbors)
    } else {
      node <- get_parent(tree, node)
    }
  }
}

#' Check if a child is a tip
#'
#' @param tree an object of class phylo
#' @param node the node to check
#' @return true or false
is_tip <- function(tree, node) {
    children <- get_children(tree, node)
    return(length(children) ==1 && node == children)
}

#' Get all the tip children of a node.
#'
#' @param tree an object of class phylo
#' @param node the node to check
#' @return the tips
get_all_children <- function(tree, node) {
    children <- get_children(tree, node)
    tips <- c()
    
    while(length(children) > 0) {
        are_tips <- c()
        for (child in children) {
          are_tips <- append(are_tips, is_tip(tree, child))
        }
        tips <- append(tips, children[are_tips])
        
        children <- children[!are_tips]
        new_children <- c()
        for (child in children) {
            new_children <- append(new_children, get_children(tree, child))
        }
        
        children <- new_children
    }
    
    return(tips)
}

#' Tree method for getting the max child clade size of a node
#' @param tree an object of class phylo
#' @param node the node to check
#' @return the size of the largest child clade
get_max_cluster_size <- function(tree, node) {
    children <- get_children(tree, node)
    num_children <- c()
    for (child in children) {
        num_children <- append(num_children, length(get_all_children(tree, child)))
    }
    return(max(num_children))
}


#' Tree method for getting the min child clade size of a node
#' @param tree an object of class phylo
#' @param node the node to check
#' @return the size of the smallest child clade
get_min_cluster_size <- function(tree, node) {
  children <- get_children(tree, node)
  num_children <- c()
  for (child in children) {
    num_children <- append(num_children, length(get_all_children(tree, child)))
  }
  return(min(num_children))
}


#' Trivial distance function for arbitrary tree clustering
#' 
#' Number of mutations along path from tip1 to LCA(tip1, tip2)
#' Ensures if on same clade, join.
#' 
#' @param tree an object of class phylo
#' @param tip1 the first leaf
#' @param tip2 the second leaf
#' @return the trivial distance between tip1, tip2
#'
#' @export
trivial_dist <- function(tree, tip1, tip2) {
  # node depths of tree
  edges <- tree$edge
  # Get the path from tip1 to root
  path1 <- c(tip1)
  root <- find_root(tree)
  parent <- tip1
  while (T) {
    parent <- edges[, 1][edges[, 2] == parent]
    path1 <- append(path1, parent)
    if (parent == root) {
      break
    }
  }
  
  mrca <- getMRCA(tree, c(tip1, tip2)) # MRCA of both
  
  
  # Depths of the internal nodes that represent the parents of tip1, tip2 right before LCA
  # ie the 'diverging' point or first split
  path_length <- which(path1 == mrca)
  
  # Return the absolute difference between the depths
  return(path_length)
}

#' Depth of tip1 parent immediately after LCA(tip1, tip2)
#' 
#' @param tree an object of class phylo
#' @param tip1 the first leaf
#' @param tip2 the second leaf
#' @return the trivial distance between tip1, tip2
#'
#' @export
lca_based_depth <- function(tree, tip1, tip2) {
  depths <- node.depth(tree)
  edges <- tree$edge
  # Get the path from tip1 to root
  path1 <- c(tip1)
  root <- find_root(tree)
  parent <- tip1
  while (T) {
    parent <- edges[, 1][edges[, 2] == parent]
    path1 <- append(path1, parent)
    if (parent == root) {
      break
    }
  }
  
  mrca <- getMRCA(tree, c(tip1, tip2)) # MRCA of both
  
  # Depths of the internal nodes that represent the parents of tip1, tip2 right before LCA
  # ie the 'diverging' point or first split
  lca_child_depth <- depths[path1[which(path1 == mrca) - 1]]
  
  # Return the absolute difference between the depths
  return(lca_child_depth)
}

#' Depth first traversal of tree.
#' 
#' @param tree a rooted tree of class `phylo`
#' @param postorder a boolean value indicating to return the list of nodes in postorder
#' @return an ordering of nodes from the traversal
depth_first_traverse_nodes <- function(tree, postorder=T) {
    root = find_root(tree)

    traversal = c()
    queue = c(root)
    while (length(queue) > 0) {
        node = queue[[1]]
        traversal = c(node, traversal)
        queue = queue[-1]
        children = get_children(tree, node)
        for (child in children) {
            if (!(child %in% traversal)) {
                queue = c(queue, child)
            }
        }
    }
    
    if (!postorder) {
        return(rev(traversal))
    }
    return(traversal)
}