#' Functions for generating projections
#'
#' This module handles the generation of lower-dimensional
#' projections from the higher-dimensional data objects.

require("fastICA")
require('Rtsne')
#require("Rtsne.multicore")
require('igraph')
require("rsvd")
require("RDRToolbox")
require("wordspace")
#require("profmem")
require("matrixStats")


#' Registers the projection methods to be used
#'
#' @param lean If FALSE, all projections applied; else a subset of essential ones are applied. Default is FALSE.
#' @return List of projection methods to be applied.
registerMethods <- function(lean=FALSE) {

  projMethods <- c()
  if (!lean) {
	projMethods <- c(projMethods, "ISOMap" = applyISOMap)
	projMethods <- c(projMethods, "ICA" = applyICA)
	#projMethods <- c(projMethods, "RBF Kernel PCA" = applyRBFPCA)
  }

  projMethods <- c(projMethods, "tSNE30" = applytSNE30)
  projMethods <- c(projMethods, "KNN" = applyKNN)
  projMethods <- c(projMethods, "tSNE10" = applytSNE10)
  projMethods <- c(projMethods, "PPT" = applySimplePPT)

  return(projMethods)
}

#' Projects data into 2 dimensions using a variety of linear and non-linear methods.
#'
#' @param expr ExpressionData object
#' @param weights weights estimated from FNR curve
#' @param filterName name of filter, to extract correct data for projections
#' @param inputProjections Precomputed projections
#' @param numCores Number of cores to use during dimensionality reductoin
#' @param lean If T, diminished number of algorithms applied, if FALSE all algorithms applied. Default is FALSE
#' @param perm_wPCA If TRUE, apply permutation wPCA to determine significant number of components. Default is FALSE.
#' @return List of Projection objects mapping projection type to projection data
#' @return List of gene names used during the dimensionality reduction phase
#' @return List of relevant PPT parameters.
generateProjections <- function(expr, weights, filterName="", inputProjections=c(), numCores = 1, lean=FALSE, perm_wPCA=FALSE) {

  if (filterName == "novar") {
    exprData <- expr@noVarFilter
  } else if (filterName == "threshold") {
    exprData <- expr@thresholdFilter
  } else if (filterName == "fano") {
    exprData <- expr@fanoFilter
  } else if (filterName == "") {
    exprData <- expr@data
  } else {
    stop("FilterName not recognized: ", filterName)
  }

  methodList = registerMethods(lean)
  t <- Sys.time()
  timingList <- (t - t)
  timingNames <- c("Start")

  #if (lean) {
  #  print("PCA")
  #  pca_res <- applyPCA(exprData, N=30)[[1]]
  #  proj <- Projection("PCA: 1,2", t(pca_res[c(1,2),]))
  #  inputProjections <- c(inputProjections, proj)
  #  inputProjections <- c(inputProjections, Projection("PCA: 1,3", t(pca_res[c(1,3),])))
  #  inputProjections <- c(inputProjections, Projection("PCA: 2,3", t(pca_res[c(2,3),])))
  #} else {
    if (perm_wPCA) {
      res <- applyPermutationWPCA(exprData, weights, components=30)
      pca_res <- res[[1]]
      print(dim(pca_res))
      loadings <- res[[3]]
    } else {
      res <- applyWeightedPCA(exprData, weights, maxComponents = 30)
      pca_res <- res[[1]]
      loadings <- res[[3]]
      #m <- profmem(pca_res <- applyWeightedPCA(exprData, weights, maxComponents = 30)[[1]])
    }
    proj <- Projection("PCA: 1,2", t(pca_res[c(1,2),]))
    inputProjections <- c(inputProjections, proj)
    inputProjections <- c(inputProjections, Projection("PCA: 1,3", t(pca_res[c(1,3),])))
    inputProjections <- c(inputProjections, Projection("PCA: 2,3", t(pca_res[c(2,3),])))
    fullPCA <- t(pca_res[1:min(15,nrow(pca_res)),])
    pca_res <- pca_res[1:5,]

    fullPCA <- as.matrix(apply(fullPCA, 2, function(x) return( x - mean(x) )))

    r <- apply(fullPCA, 1, function(x) sum(x^2))^(0.5)
    r90 <- quantile(r, c(.9))[[1]]

    if (r90 > 0) {
      fullPCA <- fullPCA / r90
    }
    fullPCA <- t(fullPCA)

  timingList <- rbind(timingList, c(difftime(Sys.time(), t, units="sec")))
  timingNames <- c(timingNames, "wPCA")
  #memProf <- c(total(m)/1e6)
  #memNames <- c("wPCA")

  PPT <- list()

  for (method in names(methodList)){
  	gc()
    message(method)
    #m <- profmem({
    if (method == "ICA" || method == "RBF Kernel PCA") {
      res <- methodList[[method]](exprData)
      proj <- Projection(method, res)
      inputProjections <- c(inputProjections, proj)
    } else {
      res <- methodList[[method]](pca_res, numCores)
      if (method == "PPT") {
      	c <- res[[1]]
      	ncls <- apply(res[[3]], 1, which.min)
        PPT <- list(res[[1]], res[[2]], res[[3]], ncls)
		PPT_neighborhood <- findNeighbors(pca_res, PPT[[1]], ncol(pca_res) / ncol(PPT[[1]]),  numCores)
      } else {
        proj <- Projection(method, res)
        inputProjections <- c(inputProjections, proj)
      }
    }
    #})
    #memProf <- rbind(memProf, total(m)/1e6)
    #memNames <- c(memNames, method)

	timingList <- rbind(timingList, c(difftime(Sys.time(), t, units="sec")))
	timingNames <- c(timingNames, method)
  }
  rownames(timingList) <- timingNames
  #return(timingList)
  #rownames(memProf) <- memNames
  #return(memProf)

  output <- list()

  for (p in inputProjections) {
    coordinates <- p@pData

    coordinates <- as.matrix(apply(coordinates, 2, function(x) return( x - mean(x) )))

    r <- apply(coordinates, 1, function(x) sum(x^2))^(0.5)
    r90 <- quantile(r, c(.9))[[1]]

    if (r90 > 0) {
      coordinates <- coordinates / r90
    }

    coordinates <- t(coordinates)
    p <- updateProjection(p, data=coordinates)
    output[[p@name]] = p

  }

  # Readjust coordinates of PPT
  c <- t(PPT[[1]][c(1,2),])
  coord <- as.matrix(apply(c, 2, function(x) return(x - mean(x))))
  r <- apply(coord, 1, function(x) sum(x^2))^(0.5)
  r90 <- quantile(r, c(0.9))[[1]]
  if (r90 > 0) {
  	  coord <- coord / r90
  }
  coord <- t(coord)
  PPT[[1]] <- coord

  output2 <- list()
  # Reposition tree node coordinatates in nondimensional space
  for (p in output) {
	if (p@name == "tSNE30" || p@name == "tSNE10" || p@name == "ISOMap") {
		new_coords <- matrix(sapply(PPT_neighborhood, function(n)  {
			n_vals <- p@pData[,n]
			centroid <- apply(n_vals, 1, mean)
			return(centroid)
		}), nrow=2, ncol=length(PPT_neighborhood))
		p@PPT_C <- new_coords
	} else {
		p@PPT_C <- PPT[[1]]
	}

	output2[[p@name]] = p
  }

  output2[["KNN"]]@PPT_C <- output2[["tSNE30"]]@PPT_C

  return(list(output2, rownames(exprData), PPT[[2]], fullPCA, loadings))
}

#' Performs weighted PCA on data
#'
#' @param exprData Expression matrix
#' @param weights Weights to use for each coordinate in data
#' @param maxComponents Maximum number of components to calculate
#' @return Weighted PCA data
#' @return Variance of each component
#' @return Eigenvectors of weighted covariance matrix, aka the variable loadings
applyWeightedPCA <- function(exprData, weights, maxComponents=200) {

  message("Weighted PCA")
  set.seed(RANDOM_SEED)

  projData <- exprData
  if (nrow(projData) != nrow(weights) || ncol(projData) != ncol(weights)) {
    weights <- weights[rownames(exprData), ]
  }

  # Center data
  wmean <- as.matrix(rowSums(multMat(projData, weights)) / rowSums(weights))
  dataCentered <- as.matrix(apply(projData, 2, function(x) x - wmean))

  # Compute weighted data
  wDataCentered <- multMat(dataCentered, weights)
  #wDataCentered <- dataCentered * weights

  # Weighted covariance / correlation matrices
  W <- tcrossprod(wDataCentered)
  Z <- tcrossprod(weights)

  wcov <- W / Z
  wcov[which(is.na(wcov))] <- 0.0
  var <- diag(wcov)

  # SVD of wieghted correlation matrix
  ncomp <- min(ncol(projData), nrow(projData), maxComponents)
  decomp <- rsvd::rsvd(wcov, k=ncomp)
  evec <- t(decomp$u)

  # Project down using computed eigenvectors
  dataCentered <- dataCentered / sqrt(var)
  wpcaData <- crossprod(t(evec), dataCentered)
  wpcaData <- wpcaData * (decomp$d*decomp$d)
  eval <- as.matrix(apply(wpcaData, 1, var))
  totalVar <- sum(apply(projData, 1, var))
  eval <- eval / totalVar

  colnames(evec) <- rownames(exprData)


  return(list(wpcaData, eval, t(evec)))

}

#' Applies pemutation method to return the most significant components of weighted PCA data
#'
#' @param expr Expression data
#' @param weights Weights to apply to each coordinate in data
#' @param components Maximum components to calculate. Default is 50.
#' @param p_threshold P Value to cutoff components at. Default is .05.
#' @param verbose Logical value indicating whether or not a verbose session is being run.
#' @return Weighted PCA data
#' @return Variance of each component.
#' @return Eigenvectors of the weighted covariance matrix
applyPermutationWPCA <- function(expr, weights, components=50, p_threshold=.05, verbose=FALSE) {
  #' Computes weighted PCA on data. Returns only significant components.
  #'
  #' After performing PCA on the data matrix, this method then uses a permutation
  #' procedure based on Buja A and Eyuboglu N (1992) to asses components for significance.
  #'
  #' Parameters:
  #'  data: (Num_Features x Num_Samples) matrix
  #'    Matrix containing data to project
  #'  weights: (Num_Features x Num_Samples) matrix
  #'    Matrix containing weights to use for each coordinate in data
  #'  components: numerical
  #'    Max components to calculate
  #'  p_threshold: numerical
  #'    P-value to cutoff components at
  #'  verbose: logical
  #'
  #'  Return:
  #		Weighed PCA data
  #		Variance of each component
  #		Eigenvectors of weighted covariance matrix.

  message("Permutation WPCA")
  comp <- min(components, nrow(expr), ncol(data))

  NUM_REPEATS <- 20;

  w <- applyWeightedPCA(expr, weights, comp)
  wPCA <- w[[1]]
  eval <- w[[2]]
  evec <- w[[3]]

  # Instantiate matrices for background distribution
  bg_vals <- matrix(0L, nrow=NUM_REPEATS, ncol=components)
  bg_data <- matrix(0L, nrow=nrow(expr), ncol=ncol(expr))
  bg_weights <- matrix(0L, nrow=nrow(expr), ncol=ncol(expr))

  # Compute background data and PCAs for comparing p values
  for (i in 1:NUM_REPEATS) {
    for (j in 1:nrow(expr)) {
      random_i <- sample(ncol(expr));
      bg_data[j,] <- expr[j,random_i]
      bg_weights[j,] <- weights[j,random_i]
    }

    bg = applyWeightedPCA(bg_data, bg_weights, comp)
    bg_vals[i,] = bg[[2]]
  }

  mu <- as.matrix(apply(bg_vals, 2, mean))
  sigma <- as.matrix(apply(bg_vals, 2, biasedVectorSD))
  sigma[which(sigma==0)] <- 1.0

  # Compute pvals from survival function & threshold components
  pvals <- 1 - pnorm((eval - mu) / sigma)
  thresholdComponent_i = which(pvals > p_threshold, arr.ind=TRUE)
  if (length(thresholdComponent_i) == 0) {
    thresholdComponent <- nrow(wPCA)
  } else {
    thresholdComponent <- thresholdComponent_i[[1]]
  }

  if (thresholdComponent < 5) {
    message("Less than 5 components identified as significant.  Preserving top 5.")
    thresholdComponent <- 5
  }

  wPCA <- wPCA[1:thresholdComponent, ]
  eval <- eval[1:thresholdComponent]
  evec = evec[1:thresholdComponent, ]

  return(list(wPCA, eval, evec))
}

#' Performs PCA on data
#'
#' @param exprData Expression data
#' @param N Number of components to retain. Default is 0
#' @param variance_proportion Retain top X PC's such that this much variance is retained; if N=0, then apply this method
#' @return Matrix containing N components for each sample.
applyPCA <- function(exprData, N=0, variance_proportion=1.0) {
  #' Performs PCA on data
  #'
  #' Args:
  #'  data: (Num_Features x Num_Samples) matrix
  #'    Matrix containing data to project into 2D
  #'  N: int
  #'    Number of Principle Components to reatin
  #'  variance_proportion: float
  #'    Retain top X principal components such taht a total of <variance_proportion> of the
  #'    variance is retained
  #'
  #' Returns:
  #'  pca_data: (Num_Components x Num_Samples) matrix
  #'    Data transformed using PCA. Num_Components = Num_Samples


  dataT = t(exprData)

  res <- prcomp(x=dataT, retx=TRUE, center=TRUE)

  if(N == 0) {

    total_var <- as.matrix(cumsum(res$sdev^2 / sum(res$sdev^2)))
    last_i <- tail(which(total_var <= variance_proportion), n=1)
    N <- last_i
  }
  return(list(t(res$x[,1:N])*-1, t(res$rotation)))
}

#' Performs ICA on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param numCores Number of cores to use
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyICA <- function(exprData, numCores) {



  set.seed(RANDOM_SEED)


  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- fastICA(ndataT, n.comp=2, maxit=100, tol=.00001, alg.typ="parallel", fun="logcosh", alpha=1,
                 method = "C", row.norm=FALSE, verbose=TRUE)

  res <- res$S
  rownames(res) <- colnames(exprData)

  return(res)
}

#' Performs Spectral Embedding  on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param numCores Number of cores to use
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applySpectralEmbedding <- function(exprData, numCores) {

  set.seed(RANDOM_SEED)

  ndata <- colNormalization(exprData)

  adj <- as.matrix(dist.matrix(t(ndata)))
  adm <- graph_from_adjacency_matrix(adj, weighted=T)
  res <- embed_adjacency_matrix(adm, 2)$X

  rownames(res) <- colnames(exprData)

  return(res)

}

#' Performs tSNE with perplexity 10 on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param numCores Number of cores to use
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applytSNE10 <- function(exprData, numCores) {

  set.seed(RANDOM_SEED)

  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  #res <- Rtsne.multicore(ndataT, dims=2, max_iter=600, perplexity=10.0, check_duplicates=F, pca=F, num_threads=numCores)
  res <- Rtsne(ndataT, dims=2, max_iter=800, perplexity=10.0, check_duplicates=F, pca=F)
  res <- res$Y
  rownames(res) <- colnames(exprData)
  return(res)

}

#' Performs tSNE with perplexity 30 on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param numCores Number of cores to use
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applytSNE30 <- function(exprData, numCores) {

  set.seed(RANDOM_SEED)

  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  #res <- Rtsne.multicore(ndataT, dims=2, max_iter=600,  perplexity=30.0, check_duplicates=F, pca=F, num_threads=numCores)
  res <- Rtsne(ndataT, dims=2, max_iter=800, perplexity=30.0, check_duplicates=F, pca=F)
  res <- res$Y

  rownames(res) <- colnames(exprData)

  return(res)
}


applyKNN <- function(exprData, numCores) {

	set.seed(RANDOM_SEED)

	k <- ball_tree_knn(t(exprData), 30, numCores)
	nn <- k[[1]]
	d <- k[[2]]

	sigma <- apply(d, 1, max)

    sparse_weights <- exp(-1 * (d * d) / sigma^2)
    weights <- load_in_knn(nn, sparse_weights)

    weightsNormFactor <- Matrix::rowSums(weights)
    weightsNormFactor[weightsNormFactor == 0] <- 1.0
    weightsNormFactor[is.na(weightsNormFactor)] <- 1.0
    weights <- weights / weightsNormFactor

    return(weights)

}

#' Performs ISOMap on data
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param numCores Number of cores to use
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyISOMap <- function(exprData, numCores) {

  set.seed(RANDOM_SEED)

  res <- Isomap(t(exprData), dims=2)
  res <- res$dim2

  rownames(res) <- colnames(exprData)

  return(res)

}

#' Performs PCA on data that has been transformed with the Radial Basis Function.
#'
#' @param exprData Expression data, NUM_GENES x NUM_SAMPLES
#' @param numCores Number of cores to use
#' @return Reduced data NUM_SAMPLES x NUM_COMPONENTS
applyRBFPCA <- function(exprData, numCores) {

  set.seed(RANDOM_SEED)

  ndata <- colNormalization(exprData)
  distanceMatrix <- as.matrix(dist.matrix(t(ndata)))
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

#' Applies the Simple PPT algorithm onto the expression data.
#'
#' @param exprData Expression data -- Num_Genes x Num_Samples
#' @param numCores Number of cores to use during this analysis
#' @param nNodes Number of nodes to find. Default is sqrt(N)
#' @param sigma regularization parameter for soft-assignment of data points to nodes, used as the variance
#'        of a guassian kernel. If 0, this is estimated automatically
#' @param gamma graph-level regularization parameter, controlling the tradeoff between the noise-levels
#'        in the data and the graph smoothness. If 0, this is estimated automatically.
#' @return Positions of the nodes in NUM_FEATURES dimensional space
#' @return Binary adjacency matrix of fitted tree
#' @return Distance matrix between fitted tree nodes and data points
#' @return Score indicating how strucutred the data is, as the z-score of the data MSE from shuffled MSE

applySimplePPT <- function(exprData, numCores, nNodes_ = round(sqrt(ncol(exprData))), sigma=0, gamma=0) {

  #exprData <- t(exprData)

  MIN_GAMMA <- 1e-5
  MAX_GAMMA <- 1e5
  DEF_TOL <- 1e-2
  DEF_MAX_ITER <- 50

  C <- NULL
  Wt <- NULL

  if (sigma == 0) {
    km <- kmeans(t(exprData), centers=round(sqrt(ncol(exprData))), nstart=1, iter.max=50)$centers

    sigma <- mean(apply(as.matrix(sqdist(t(exprData), km)), 1, min))
  }

  if (gamma == 0) {

    currGamma <- MIN_GAMMA
    nNodes <- round(log(ncol(exprData)))

    prevMSE <- -Inf
    minMSE <- Inf
    minMSEGamma <- MIN_GAMMA

    tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
    C <- tr[[1]]
    Wt <- tr[[2]]
    currMSE <- tr[[3]]

    while ( ((prevMSE / currMSE) - 1 < 0.05) && currGamma <= MAX_GAMMA) {
      prevMSE <- currMSE
      currGamma <- currGamma * 10
      tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
      C <- tr[[1]]
      Wt <- tr[[2]]
      currMSE <- tr[[3]]
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
      C <- tr[[1]]
      Wt <- tr[[2]]
      currMSE <- tr[[3]]
    }

    if (nNodes_ > nNodes) {

      if (nNodes_ != 0) {
        nNodes <- nNodes_
      } else {
        nNodes <- round(sqrt(ncol(exprData)))
      }

      tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
      C <- tr[[1]]
      Wt <- tr[[2]]
      currMSE <- tr[[3]]

      # Calculate the degree distribution
      deg <- colSums(Wt)
      br <- seq(0, max(1, max(deg)))
      degDist <- hist(deg, br, plot=F)$counts
      deg_g2c <- 0
      if (length(degDist) > 2) {
        deg_g2c <- sum(degDist[seq(3, length(degDist))])
      }
      deg_g2f <- deg_g2c / nNodes

      while ( !(deg_g2c > 0 && (deg_g2f <= 0.1 || deg_g2c < 5)) && (currGamma >= minGamma) ) {
        currGamma <- currGamma / sqrt(10)
        tr <- fitTree(exprData, nNodes, sigma, currGamma, DEF_TOL, DEF_MAX_ITER)
        C <- tr[[1]]
        Wt <- tr[[2]]
        currMSE <- tr[[3]]

        deg <- colSums(Wt)
        br <- seq(0, max(1, max(deg)))
        degDist <- hist(deg, br, plot=F)$counts
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
  C <- tr[[1]]
  Wt <- tr[[2]]
  mse <- tr[[3]]


  return(list(C, Wt, sqdist(t(exprData), t(C)), mse))
}

#' Fit tree using input parameters
#'
#' @param expr Data to fit (NUM_GENES x NUM_SAMPLES)
#' @param nNodes Number of nodes in the fitted tree, default is square-root of number of data points
#' @param sigma Regularization parameter for soft-assignment of data points to nodes, used as the
#'              variance of a gaussian kernel. If 0, this is estimated automatically.
#' @param gamma Graph-level regularization parameter, controlling the tradeoff between the noise-levels
#'              in the data and the graph smoothness. If 0, this is estimated automatically
#' @param tol Tolerance to use when fitting the tree
#' @param maxIter Maximum number of Iterations ot run the algorithm for
#' @return Positions of the nodes in NUM_FEATURES dimensional space
#' @return Binary adjacency matrix of fitted tree
#' @return MSE between Fitted nodes and original data points

fitTree <- function(expr, nNodes, sigma, gamma, tol, maxIter) {
  #' Fit tree using input parameters
  #'
  #' Paramters:
  #'  expr: (NUM_GENES x NUM_SAMPLES) matrix
  #'    data to fit
  #'  nNodes: numeric
  #'    number of nodes in the fitted tree, default is the square-root of number of data points
  #'  sigma: numeric
  #'    regularization parameter for soft-assignment of data points to nodes, used as the
  #'    variance of a guassian kernel. If 0, this is estimated automatically
  #'  gamma: numeric
  #'    graph-level regularization parameter, controlling the tradeoff between the noise-levels
  #'    in the data and the graph smoothness. If 0, this is estimated automatically.


  km <- kmeans(t(expr), centers=nNodes, nstart=10, iter.max=100)$centers
  cc_dist <- as.matrix(sqdist(km, km))
  cx_dist <- as.matrix(sqdist(t(expr), km))
  prevScore = Inf
  currScore = -Inf
  currIter = 0

  while (!(prevScore - currScore < tol) && !(currIter > maxIter)) {
    currIter <- currIter + 1
    prevScore <- currScore
    W <- mst(graph_from_adjacency_matrix(cc_dist, weighted= T, mode="undirected"))
    Wt <- get.adjacency(W, sparse=FALSE)

    Ptmp <- -(cx_dist / sigma)
    Psums <- matrix(rep(apply(Ptmp, 1, logSumExp), each=ncol(Ptmp)), nrow=nrow(Ptmp), ncol=ncol(Ptmp), byrow=T)
    P <- exp(Ptmp - Psums)

    delta <- diag(colSums(P))
    L <- laplacian_matrix(W)
    xp <- crossprod(t(expr), P)
    invg <- as.matrix(solve( ((2 / gamma) * L) + delta))
    C <- tcrossprod(xp, invg)

    cc_dist <- as.matrix(dist.matrix(t(C)))
    cx_dist <- as.matrix(sqdist(t(expr), t(C)))

    P <- clipBottom(P, mi=min(P[P>0]))
    currScore <- sum(Wt * cc_dist) + (gamma * sum(P * ((cx_dist) + (sigma * log(P)))))

  }

  return(list(C, Wt, getMSE(C, expr)))

}


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
	x[which(x<0)] <- 0
	return(x)

}

#' Sets all values below a certain level in the data equal to 0
#'
#' @param x Data matrix
#' @param mi Minimum value
#' @return Data matrix with all values less than MI set to 0
clipBottom <- function(x, mi) {
  # Sets all values below MI to 0.
  #
  # Args:
  #		x: data matrix
  #		mi: mininum value
  #
  # Returns:
  #		Data matrix with all values less than MI set to 0
  x[x < mi] <- mi
  return(x)
}

findNeighbors <- function(data, query, k, numCores) {

	neighborhood <- lapply(1:ncol(query), function(x) {
		vkn <- ball_tree_vector_knn(t(data), query[,x], k, numCores)
		return(vkn[[1]])
	})

	return(neighborhood)
}

