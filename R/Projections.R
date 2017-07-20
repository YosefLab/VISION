#' Functions for generating projections
#' 
#' This module handles the generation of lower-dimensional
#' projections from the higher-dimensional data objects.

require("fastICA")
require('Rtsne')
require("Rtsne.multicore")
require('igraph')
require("rsvd")
require("RDRToolbox")
require("wordspace")
#require("profmem")


registerMethods <- function(lean=FALSE) {
  
  projMethods <- c()
  if (!lean) {
	projMethods <- c(projMethods, "ISOMap" = applyISOMap)
	projMethods <- c(projMethods, "ICA" = applyICA)
	#projMethods <- c(projMethods, "RBF Kernel PCA" = applyRBFPCA)
  }
  
  projMethods <- c(projMethods, "tSNE30" = applytSNE30)
  projMethods <- c(projMethods, "tSNE10" = applytSNE10)
  #projMethods <- c(projMethods, "PPT" = applySimplePPT)
  
  return(projMethods)
}

generateProjections <- function(expr, weights, filterName="", inputProjections=c(), numCores = 1, lean=FALSE, perm_wPCA=FALSE, optClust=0, approximate=F) {
  #' Projects data into 2 dimensions using a variety of linear and non-linear methods
  #' 
  #' Parameters:
  #'  expr: (ExpressionData) expression matrix to project into 2D
  #'  filterName: (character) name of filter to apply to signatures, should match a filter that was
  #'              applied to the exprData before
  #'  inputProjections: (list of ProjectionData) collection of projection data types; names are of type
  #'              character, mapping to a projection             
  #' 
  #' Returns:
  #'  projections: (list of Projection) collection mapping projection type to the actual projection
  #'  PC_data: (matrix) weighted PCA of original data type
  
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
  #  pca_res <- applyPCA(exprData, N=30)
  #  proj <- Projection("PCA: 1,2", t(pca_res[c(1,2),]))
  #  inputProjections <- c(inputProjections, proj)
  #  inputProjections <- c(inputProjections, Projection("PCA: 1,3", t(pca_res[c(1,3),])))
  #  inputProjections <- c(inputProjections, Projection("PCA: 2,3", t(pca_res[c(2,3),])))
  #} else {
    if (perm_wPCA) {
      pca_res <- applyPermutationWPCA(exprData, weights, components=30)[[1]]
    } else {
      pca_res <- applyWeightedPCA(exprData, weights, maxComponents = 30)[[1]]
      #m <- profmem(pca_res <- applyWeightedPCA(exprData, weights, maxComponents = 30)[[1]])
	}
    proj <- Projection("PCA: 1,2", t(pca_res[c(1,2),]))
    inputProjections <- c(inputProjections, proj)
    inputProjections <- c(inputProjections, Projection("PCA: 1,3", t(pca_res[c(1,3),])))
    inputProjections <- c(inputProjections, Projection("PCA: 2,3", t(pca_res[c(2,3),])))
  #}
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
      	ncls <- apply(res[[3]], 1, which.min)
        PPT <- list(res[[1]], res[[2]], res[[3]], ncls)
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
  
  return(list(output, rownames(exprData), PPT))
}

applyWeightedPCA <- function(exprData, weights, maxComponents=200) {
  #' Performs Weighted PCA on the data
  #' 
  #' Parameters:
  #'  data: (Num_Features x Num_Samples) matrix
  #'    Matrix containing data to project 
  #'  weights: (Num_Features x Num_Samples) matrix
  #'    Matrix containing weights to use for each coordinate in data
  #'  max_components: numeric
  #'    Maximum number of components to calculate
  #'    
  #' Returns:
  #'  pca_data: (Num_Components x Num_Samples) matrix
  #'    Data transformed using PCA.  Num_Components = Num_Samples
  
  message("Weighted PCA")
  set.seed(RANDOM_SEED)
  
  projData <- exprData
  if (nrow(projData) != nrow(weights) || ncol(projData) != ncol(weights)) {
    weights <- weights[rownames(exprData), ]
  }
  
  # Center data
  wmean <- as.matrix(rowSums(multMat(projData, weights)) / rowSums(weights))
  #wmean <- as.matrix(rowSums(projData * weights) / rowSums(weights))
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
  # NOTE: Weighted Covariance works better than Weighted Correlation for computing the eigenvectors
  # NOTE: rsvd() method is 5x more memory efficent than eigs(); eigs() is slightly faster  
  decomp <- rsvd::rsvd(wcov, k=ncomp)
  evec <- t(decomp$u)
  
  # Project down using computed eigenvectors
  dataCentered <- dataCentered / sqrt(var)
  wpcaData <- crossprod(t(evec), dataCentered)
  eval <- as.matrix(apply(wpcaData, 1, var))
  totalVar <- sum(apply(projData, 1, var))
  eval <- eval / totalVar

  
  return(list(wpcaData, eval, t(evec)))
    
}

applyPermutationWPCA <- function(expr, weights, components=50, p_threshold=.05, verbose=FALSE, debug=FALSE) {
  #' Computes weighted PCA on data. Returns only significant components.
  #' 
  #' After performing PCA on the data matrix, this method then uses a permutation
  #' procedure based on Buja A and Eyuboglu N (1992) to asses components for significance.
  #' 
  #' Paramaters:
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
  #'    reducedData: (Num_Components X Num_Samples)
  
  message("Permutation WPCA")
  comp <- min(components, nrow(expr), ncol(data))
  
  NUM_REPEATS <- 1;
  
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
    
    print(i)
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
applyPCA <- function(exprData, N=0, variance_proportion=1.0) {
  #' Performs PCA on data
  #' 
  #' Parameters:
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
  return(t(res$x[,1:N])*-1)
}

applyICA <- function(exprData, numCores, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- fastICA(ndataT, n.comp=2, maxit=100, tol=.00001, alg.typ="parallel", fun="logcosh", alpha=1,
                 method = "C", row.norm=FALSE, verbose=TRUE)
  
  return(res$S)
}

applySpectralEmbedding <- function(exprData, numCores, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  
  adj <- as.matrix(dist.matrix(t(ndata)))
  adm <- graph_from_adjacency_matrix(adj, weighted=T)
  res <- embed_adjacency_matrix(adm, 2)$X
  
  rownames(res) <- colnames(exprData)
 
  return(res)
  
}

applytSNE10 <- function(exprData, numCores, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- Rtsne.multicore(ndataT, dims=2, max_iter=600, perplexity=10.0, check_duplicates=F, pca=F, num_threads=numCores)
  res <- res$Y
  rownames(res) <- colnames(exprData)
  return(res)
  
}

applytSNE30 <- function(exprData, numCores, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- Rtsne.multicore(ndataT, dims=2, max_iter=600,  perplexity=30.0, check_duplicates=F, pca=F, num_threads=numCores)
  res <- res$Y
  
  rownames(res) <- colnames(exprData)
  
  return(res)
}

applyISOMap <- function(exprData, numCores, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  res <- Isomap(t(exprData), dims=2)
  res <- res$dim2
  
  rownames(res) <- colnames(exprData)
  
  return(res)
  
}

applyRBFPCA <- function(exprData, numCores, projWeights=NULL) {
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


applySimplePPT <- function(exprData, numCores, projWeights=NULL, nNodes_ = 0, sigma=0, gamma=0) {
  #' Principle Tree Analysis
  #' 
  #' After begin initialized by either the fit or autoFit functions, the resulting tree structure
  #' is contained within these fields:
  #'  C: (NUM_FEATURES x NUM_NODES) matrix
  #'    positions of the nodes in NUM_FEATURES dimensional space
  #'  W: (NUM_NODES X NUM_NODES) matrix
  #'    binary adjacency matrix of the fitted tree
  #'  Distmat: (NUM_NODES X NUM_SAMPLES) matrix
  #'    distance matrix between fitted tree nodes and the data points
  #'  structScore: (numeric)
  #'    score indicating how structured the data is, as the z-score of the data MSE from shuffled MSE
  
  exprData <- t(exprData)
  
  MIN_GAMMA <- 1e-5
  MAX_GAMMA <- 1e5
  DEF_TOL <- 1e-3
  DEF_MAX_ITER <- 50
  
  C <- NULL
  Wt <- NULL
  
  if (nNodes_ == 0) {
    nNodes_ <- round(sqrt(ncol(exprData)))
  }
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
    minGamma <- MIN_GAMMA
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

fitTree <- function(expr, nNodes, sigma, gamma, tol, maxIter) {
  #' Fit tree using input parameters
  #' 
  #' Paramters:
  #'  expr: (NUM_GENES x NUM_SAMPLES) matrix
  #'    data to fit
  #'  nNodes: numeric
  #'    number of nodes in the fitted tree, default is the square-root of number of data points
  #'  sigma: numeric
  #'    regularization paramter for soft-assignment of data points to nodes, used as the
  #'    variance of a guassian kernel. If 0, this is estimated automatically
  #'  gamma: numeric
  #'    graph-level regularization parameter, controlling the tradeoff between the noise-levels
  #'    in the data and the graph smoothness. If 0, the is estimated automatically.
  
  
  km <- kmeans(t(expr), centers=nNodes, nstart=10, iter.max=100)$centers
  cc_dist <- as.matrix(dist.matrix(km))
  cx_dist <- as.matrix(sqdist(t(expr), km))
  prevScore = Inf
  currScore = 0
  currIter = 0
  
  while (!( (prevScore - currScore < tol) || (currIter > maxIter) )){
    currIter <- currIter + 1
    prevScore <- currScore
    W <- mst(graph_from_adjacency_matrix(cc_dist, weighted=T, mode="undirected"))
    Wt <- get.adjacency(W, sparse=FALSE)
    
    Ptmp <- -(cx_dist / sigma)
    Psums <- matrix(rep(apply(Ptmp, 1, logSumExp), each=ncol(Ptmp)), ncol=ncol(Ptmp), byrow=T)
    P <- exp(Ptmp - Psums)
    
    delta <- diag(colSums(P))
    L <- laplacian_matrix(W)
    xp <- crossprod(t(expr), P)
    invg <- as.matrix(solve( ((2 / gamma) * L) + delta))
    C <- crossprod(t(xp), invg)
    
    cc_dist <- as.matrix(dist.matrix(t(C)))
    cx_dist <- as.matrix(sqdist(t(expr), t(C)))
    
    P <- clipBottom(P, mi=min(P[P>0]))
    currScore <- sum(Wt * cc_dist) + (gamma * sum(P * ((cx_dist) + (sigma * log(P)))))
    
  }
  
  return(list(C, Wt, getMSE(C, expr)))
  
}

getMSE <- function(C, X) {
  if (is.na(C) || is.na(X)) {
    return(NULL)
  }
  mse <- mean( apply( as.matrix(sqdist(t(X), t(C))), 1, min))
  return(mse)
  
}

sqdist <- function(X, Y) {
  #' Alternative computation of distance matrix, based on matrix multiplication.
  #' X is n x d matrix
  #' Y is m x d matrix
  #' Returns n x m distance matrix

	aa = rowSums(X**2)
	bb = rowSums(Y**2)
	x = -2 * tcrossprod(X, Y)
	x = x + aa
	x = t(t(x) + bb)
	x[which(x<0)] <- 0
	return(sqrt(x))

}

clipBottom <- function(x, mi) {
  x[x < mi] <- mi
  return(x)
}

