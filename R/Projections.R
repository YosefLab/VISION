
#' Registers the projection methods to be used
#'
#' @param numCells number of cells in this analysis
#' @param projection_methods character vector of projection methods to use
#' @return List of projection methods to be applied.
registerMethods <- function(projection_methods, numCells) {

    all_proj_methods <- c("ISOMap" = applyISOMap,
                        "ICA" = applyICA,
                        "tSNE30" = applytSNE30,
                        "tSNE10" = applytSNE10,
                        "RBFPCA" = applyRBFPCA,
                        "UMAP" = applyUMAP)

    projMethods <- all_proj_methods[projection_methods]

    return(projMethods)
}

#' Projects data into 2 dimensions using a variety of linear and non-linear methods.
#'
#' @importFrom stats quantile
#' @param expr numeric matrix of gene expression
#' @param latentSpace numeric matrix cells x components
#' @param projection_genes character vector of gene names to use for projections
#' @param projection_methods character vector of projection methods to use
#' @param K Number of neighbors to use in projections.
#' @return list of Projection objects
generateProjectionsInner <- function(expr, latentSpace, projection_genes=NULL, projection_methods = NULL, K = round(sqrt(ncol(expr)))) {

    NUM_CELLS <- nrow(latentSpace)
    methodList <- registerMethods(projection_methods, NUM_CELLS)

    projections <- list()

    projections[["Latent Space"]] <- latentSpace

    N <- length(methodList)
    for (i in seq_len(N)){
        method <- names(methodList)[i]
        message(
            sprintf("  Running method %i/%i: %s ...", i, N, method)
        )
        ## run on raw data
        if (method == "ICA" || method == "RBFPCA") {

            if (!is.null(projection_genes)) {
                exprData <- expr[projection_genes, , drop = FALSE]
            } else {
                exprData <- expr
            }

            exprData <- matLog2(exprData)

            res <- methodList[[method]](exprData)
            projections[[method]] <- res
        } else if (method == "UMAP") {
            res <- methodList[[method]](t(latentSpace), K)
            projections[[method]] <- res
        } else { ## run on reduced data
            res <- methodList[[method]](t(latentSpace))
            projections[[method]] <- res
        }
    }

    return(projections)
}


#' Performs weighted PCA on data
#'
#' @importFrom rsvd rsvd
#'
#' @param exprData Expression matrix
#' @param weights Weights to use for each coordinate in data
#' @param maxComponents Maximum number of components to calculate
#' @return Weighted PCA data
#' @return Variance of each component
#' @return Eigenvectors of weighted covariance matrix, aka the variable loadings
applyWeightedPCA <- function(exprData, weights, maxComponents=200) {

    projData <- exprData
    if (nrow(projData) != nrow(weights) || ncol(projData) != ncol(weights)) {
        weights <- weights[rownames(exprData), ]
    }

    # Center data
    wmean <- rowSums(projData * weights) / rowSums(weights)
    dataCentered <- projData - wmean

    # Compute weighted data
    wDataCentered <- dataCentered * weights

    # Weighted covariance / correlation matrices
    W <- Matrix::tcrossprod(wDataCentered)
    Z <- Matrix::tcrossprod(weights)

    wcov <- as.matrix(W / Z)
    wcov[is.na(wcov)] <- 0.0
    var <- diag(wcov)

    # SVD of wieghted correlation matrix
    ncomp <- min(ncol(projData), nrow(projData), maxComponents)
    decomp <- rsvd::rsvd(wcov, k = ncomp)
    evec <- t(decomp$u)

    # Project down using computed eigenvectors
    wpcaData <- evec %*% dataCentered

    eval <- rowVars(wpcaData)
    totalVar <- sum(rowVars(projData))
    eval <- eval / totalVar

    colnames(evec) <- rownames(exprData)


    return(list(wpcaData, eval, t(evec)))

}

#' Applies pemutation method to return the most significant components of weighted PCA data
#'
#' @details Based on the method proposed by Buja and Eyuboglu (1992), PCA is performed on the data
#' then a permutation procedure is used to assess the significance of components
#'
#' @importFrom stats pnorm
#' @param expr Expression data
#' @param weights Weights to apply to each coordinate in data
#' @param components Maximum components to calculate. Default is 50.
#' @param p_threshold P Value to cutoff components at. Default is .05.
#' @return (list):
#' \itemize{
#'     \item wPCA: weighted PCA data
#'     \item eval: the proortinal variance of each component
#'     \item evec: the eigenvectors of the weighted covariance matrix
#'     \item permuteMatrices: the permuted matrices generated as the null distrbution
#' }
applyPermutationWPCA <- function(expr, weights, components=50, p_threshold=.05) {
    comp <- min(components, nrow(expr), ncol(expr))

    NUM_REPEATS <- 20;

    w <- applyWeightedPCA(expr, weights, comp)
    wPCA <- w[[1]]
    eval <- w[[2]]
    evec <- w[[3]]

    # Instantiate matrices for background distribution
    bg_vals <- matrix(0L, nrow=NUM_REPEATS, ncol=components)
    bg_data <- matrix(0L, nrow=nrow(expr), ncol=ncol(expr))
    bg_weights <- matrix(0L, nrow=nrow(expr), ncol=ncol(expr))

    permMats <- list()

    # Compute background data and PCAs for comparing p values
    for (i in 1:NUM_REPEATS) {
    for (j in 1:nrow(expr)) {
        random_i <- sample(ncol(expr));
        bg_data[j,] <- expr[j,random_i]
        bg_weights[j,] <- weights[j,random_i]
    }

    bg = applyWeightedPCA(bg_data, bg_weights, comp)
    bg_vals[i,] = bg[[2]]
    permMats[[i]] <- bg[[1]]
    }

    mu <- as.matrix(apply(bg_vals, 2, mean))
    sigma <- as.matrix(apply(bg_vals, 2, sd))
    sigma[sigma==0] <- 1.0

    # Compute pvals from survival function & threshold components
    pvals <- 1 - pnorm((eval - mu) / sigma)
    thresholdComponent_i = which(pvals > p_threshold, arr.ind=TRUE)
    if (length(thresholdComponent_i) == 0) {
        thresholdComponent <- nrow(wPCA)
    } else {
        thresholdComponent <- thresholdComponent_i[[1]]
    }

    if (thresholdComponent < 5) {
    # Less than 5 components identified as significant.  Preserving top 5.
    thresholdComponent <- 5
    }

    wPCA <- wPCA[1:thresholdComponent, ]
    eval <- eval[1:thresholdComponent]
    evec = evec[1:thresholdComponent, ]

    permMats <- lapply(permMats, function(m) {m[1:thresholdComponent, ]})

    return(list(wPCA = wPCA, eval = eval, evec = evec, permuteMatrices = permMats))
}

#' Performs PCA on data
#' @importFrom  utils tail
#' @importFrom stats prcomp
#' @param exprData Expression data
#' @param N Number of components to retain. Default is 0
#' @param variance_proportion Retain top X PC's such that this much variance is retained; if N=0, then apply this method
#' @return Matrix containing N components for each sample.
applyPCA <- function(exprData, N=0, variance_proportion=1.0) {
    res <- prcomp(x=t(exprData), retx=TRUE, center=TRUE)

    if(N == 0) {
    total_var <- as.matrix(cumsum(res$sdev^2 / sum(res$sdev^2)))
    last_i <- tail(which(total_var <= variance_proportion), n=1)
    N <- last_i
    }
    return(list(t(res$x[,1:N])*-1, t(res$rotation)))
}

#' Performs ICA on data
#'
#' @importFrom fastICA fastICA
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyICA <- function(exprData) {

    ndataT <- t(exprData)
    res <- fastICA(ndataT, n.comp=2, maxit=100, tol=.00001, alg.typ="parallel", fun="logcosh", alpha=1,
                    method = "C", row.norm=FALSE, verbose=TRUE)

    res <- res$S
    rownames(res) <- colnames(exprData)

    return(res)
}

#' Performs Spectral Embedding  on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @importFrom wordspace dist.matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph embed_adjacency_matrix
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applySpectralEmbedding <- function(exprData) {


    adj <- as.matrix(dist.matrix(t(exprData)))
    adm <- graph_from_adjacency_matrix(adj, weighted=TRUE)
    res <- embed_adjacency_matrix(adm, 2)$X

    rownames(res) <- colnames(exprData)

    return(res)

}

#' Performs UMAP on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param K Number of neighbors to use in UMAP projection.
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyUMAP <- function(exprData, K) {

	print(requireNamespace('uwot'))
    if (!requireNamespace("uwot", quietly = TRUE)){
        stop("Package \"uwot\" needed to run UMAP.  Please install it using:\n\n   devtools::install_github(\"jlmelville/uwot\")\n\n",
            call. = FALSE)
    }

	ndataT <- t(exprData)
    n_workers <- getOption("mc.cores")
    n_workers <- if (is.null(n_workers)) 2 else n_workers
	res <- uwot::umap(
        ndataT, n_neighbors = K,
        n_threads = n_workers, ret_nn = T
    )
	res <- res$embedding
	rownames(res) <- colnames(exprData)

	return(res)
}

#' Performs tSNE with perplexity 10 on data
#'
#' @importFrom Rtsne Rtsne
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#'
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applytSNE10 <- function(exprData) {

    ndataT <- t(exprData)
    res <- Rtsne(ndataT, dims=2, max_iter=800, perplexity=10.0,
                check_duplicates=FALSE, pca=FALSE)
    res <- res$Y
    rownames(res) <- colnames(exprData)
    return(res)

}

#' Performs tSNE with perplexity 30 on data
#'
#' @importFrom Rtsne Rtsne
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applytSNE30 <- function(exprData) {

    ndataT <- t(exprData)
    res <- Rtsne(ndataT, dims=2, max_iter=800, perplexity=30.0,
                check_duplicates=FALSE, pca=FALSE)
    res <- res$Y

    rownames(res) <- colnames(exprData)

    return(res)
}

#' Performs ISOMap on data
#'
#' @importFrom vegan isomap
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyISOMap <- function(exprData) {

    d.expr = dist.matrix(exprData, byrow=F)
    res <- isomap(d.expr, ndim=2, epsilon=1e6)
    res <- res$points

    rownames(res) <- colnames(exprData)

    return(res)

}

#' Performs PCA on data that has been transformed with the Radial Basis Function.
#'
#' @importFrom stats sd
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyRBFPCA <- function(exprData) {

    distanceMatrix <- as.matrix(dist.matrix(t(exprData)))
    distanceMatrix <- log(distanceMatrix)
    point_mult(distanceMatrix, distanceMatrix)
    kMat <- as.matrix(exp(-1 * (distanceMatrix) / .33^2))
    diag(kMat) <- 0
    kMatNormFactor <- rowSums(kMat)
    kMatNormFactor[kMatNormFactor == 0] <- 1.0
    kMatNormFactor[is.na(kMatNormFactor)] <- 1.0
    kMat <- kMat / kMatNormFactor

    # Compute normalized matrix & covariance matrix
    kMat <- as.matrix(kMat, 1, function(x) (x - mean(x)) / sd(x))
    W <- tcrossprod(kMat)

    decomp <- rsvd::rsvd(W, k=2)
    evec <- decomp$u

    # project down using evec
    rbfpca <- crossprod(t(kMat), evec)
    rownames(rbfpca) <- colnames(exprData)

    return(rbfpca)

}


#' Alternative computation of distance matrix, based on matrix multiplication.
#'
#' @param X n x d matrix
#' @param Y m x d matrix
#' @return n x m distance matrix
sqdist <- function(X, Y) {

    aa = rowSums(X**2)
    bb = rowSums(Y**2)
    x = -2 * tcrossprod(X, Y)
    x = x + aa
    x = t(t(x) + bb)
    x[x<0] <- 0
    return(x)

}

#' Sets all values below a certain level in the data equal to 0
#'
#' @param x Data matrix
#' @param mi Minimum value
#' @return Data matrix with all values less than MI set to 0
clipBottom <- function(x, mi) {
    x[x < mi] <- mi
    return(x)
}

#' compute for each vector the weights to apply to it's K nearest neighbors
#' @importFrom Matrix rowSums
#' @importFrom Matrix sparseMatrix
#' @importFrom matrixStats rowMaxs
#' @param object matrix to use for KNN
#' @param K Number of neighbors to consider.
#' @return a list of two items:
#'          indices: matrix, cells X neighbors
#'              Each row specifies indices of nearest neighbors
#'          weights: matrix, cells X neighbors
#'              Corresponding weights to nearest neighbors
setMethod("computeKNNWeights", signature(object = "matrix"),
    function(object, K = round(sqrt(ncol(object)))) {

        print(K)
        n_workers <- getOption("mc.cores")
        n_workers <- if (is.null(n_workers)) 2 else n_workers

        k <- ball_tree_knn(object, K, n_workers)
        nn <- k[[1]]
        d <- k[[2]]

        sigma <- rowMaxs(d)
        sigma[sigma == 0] <- 1.0  # Can happen if all neighbors at same point as cell
        sparse_weights <- exp(-1 * (d * d) / sigma ^ 2)

        # Normalize row sums = 1
        weightsNormFactor <- rowSums(sparse_weights)
        weightsNormFactor[weightsNormFactor == 0] <- 1.0
        sparse_weights <- sparse_weights / weightsNormFactor

        rownames(nn) <- rownames(object)
        rownames(d) <- rownames(object)

        return(list(indices = nn, weights = sparse_weights))
    }
)
