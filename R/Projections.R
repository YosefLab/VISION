#' Functions for generating projections
#' 
#' This module handles the generation of lower-dimensional
#' projections from the higher-dimensional data objects.

require("fastICA")
require('Rtsne')
require('dimRed')
require("RANN")
require('igraph')
require('kernlab')
require("vegan")
require("loe")
require("cluster")
require("smacof")
require("rARPACK")
require("BKPC")


registerMethods <- function(lean=FALSE) {
  
  #projMethods <- c("ICA" = applyICA)
  projMethods <- c()
  if (!lean) {
    projMethods <- c(projMethods, "Spectral Embedding" = applySpectralEmbedding)
    projMethods <- c(projMethods, "MDS" = applyMDS)
  }
  
  projMethods <- c(projMethods, "RBF Kernel PCA" = applyRBFPCA)
  #projMethods <- c(projMethods, "ISOMap" = applyISOMap)
  projMethods <- c(projMethods, "tSNE30" = applytSNE30)
  projMethods <- c(projMethods, "tSNE10" = applytSNE10)
  
  return(projMethods)
}

generateProjections <- function(expr, weights, filterName="", inputProjections=c(), lean=FALSE, perm_wPCA=FALSE) {
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
  } else {
    stop("FilterName not recognized: ", filterName)
  }
  
  methodList = registerMethods(lean)
  
  if (lean) {
    print("PCA")
    pca_res <- applyPCA(exprData, N=30)
    proj <- Projection("PCA: 1,2", t(pca_res[c(1,2),]))
    inputProjections <- c(inputProjections, proj)
    inputProjections <- c(inputProjections, Projection("PCA: 1,3", t(pca_res[c(1,3),])))
    inputProjections <- c(inputProjections, Projection("PCA: 2,3", t(pca_res[c(2,3),])))
  } else {
    if (perm_wPCA) {
      print("Permutation PCA")
      pca_res <- applyPermutationWPCA(exprData, weights, components=30)[[1]]
    } else {
      print("Weighted PCA")
      pca_res <- applyWeightedPCA(exprData, weights, maxComponents = 30)[[1]]
    }
    proj <- Projection("PCA: 1,2", t(pca_res[c(1,2),]))
    inputProjections <- c(inputProjections, proj)
    inputProjections <- c(inputProjections, Projection("PCA: 1,3", t(pca_res[c(1,3),])))
    inputProjections <- c(inputProjections, Projection("PCA: 2,3", t(pca_res[c(2,3),])))
  }
  
  for (method in names(methodList)){
    print(method)
    if (method == "ICA" || method == "RBF Kernel PCA") {
      res <- methodList[[method]](exprData)
      proj <- Projection(method, res)
      inputProjections <- c(inputProjections, proj)
    } else {
      # Check if method is PCA because function call isn't general
      res <- methodList[[method]](pca_res)
      proj <- Projection(method, res)
      inputProjections <- c(inputProjections, proj)
    }
  }  
  
  output <- c()
  
  for (p in inputProjections) {

    coordinates <- p@pData
    
    
    for (i in 1:ncol(coordinates)) {
      coordinates[,i] <- coordinates[,i] - mean(coordinates[,i])
    }
    
    r <- apply(coordinates, 1, function(x) sum(x^2))^(0.5)
    r90 <- quantile(r, c(.9))[[1]]

    if (r90 > 0) {
      coordinates <- coordinates / r90
    }

    coordinates <- t(coordinates)
    p <- updateProjection(p, data=coordinates)
    output <- c(output, p)
  }
  
  
  return(list(output, rownames(exprData)))
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

  set.seed(RANDOM_SEED)
        
  projData <- exprData
  if (dim(projData) != dim(weights)) {
    weights <- weights[rownames(exprData), ]
  }
  
  # Center data
  wmean <- as.matrix(apply(projData * weights, 1, sum) / apply(weights, 1, sum))
  dataCentered <- as.matrix(apply(projData, 2, function(x) x - wmean))

  # Compute weighted data
  wDataCentered <- dataCentered * weights
  
  print("wcov")
  # Weighted covariance / correlation matrices
  W <- wDataCentered %*% t(wDataCentered)
  Z <- weights %*% t(weights)
  wcov <- W / Z
  wcov[which(is.na(wcov))] <- 0.0
  var <- diag(wcov)
  cor <- wcov / sqrt(var %*% t(var))
  
  # SVD of wieghted correlation matrix
  ncomp <- min(ncol(projData), nrow(projData), maxComponents)
  print("eig")
  # NOTE: Weighted Covariance works better than Weighted Correlation for computing the eigenvectors
  eig_obj = eigs(wcov,k = ncomp,which = "LM")
  evec <- t(eig_obj$vectors)
  
  
  print("eval")
  # Project down using computed eigenvectors
  dataCentered <- dataCentered / sqrt(var)
  wpcaData <- as.matrix(evec %*% dataCentered)
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

applyICA <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- fastICA(ndataT, n.comp=2, maxit=100, tol=.00001, alg.typ="parallel", fun="logcosh", alpha=1,
                 method = "C", row.norm=FALSE, verbose=TRUE)
  
  return(res$S)
}

applySpectralEmbedding <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  
  #adj <- make.distmat(ndataT)
  adj <- as.matrix(dist(ndataT, method="euclidean"))
  res <- spec.emb(adj, 2, norm=TRUE)
  #adm <- graph_from_adjacency_matrix(adj)

  #res <- embed_adjacency_matrix(adm, 2)
  rownames(res) <- colnames(exprData)
 
  return(res)
  
}

applyMDS <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  distN <- dist(ndataT, method="euclidean")
  
  #res <- cmdscale(distN, k=2, eig=TRUE)
  res <- mds(distN, ndim=2)
  res <- res$conf
  rownames(res) <- colnames(exprData)
  return(res)
  
}

applytSNE10 <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- Rtsne(ndataT, dims=2, perplexity=10.0, pca=FALSE, theta=0.0)
  res <- res$Y
  rownames(res) <- colnames(exprData)
  return(res)
  
}

applytSNE30 <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- Rtsne(ndataT, dims=2, perplexity=30.0, pca=FALSE, theta=0.0)
  res <- res$Y
  
  rownames(res) <- colnames(exprData)
  
  return(res)
}

applyISOMap <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  d <- dist(t(ndata), method="euclidean")
  res <- isomap(d, k=4, ndim=2, fragmentedOK=TRUE)
  res <- res$points
  
  return(res)
  
}

applyRBFPCA <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- kpca(ndataT, features=2, kpar=list(sigma=0.2), kernel='rbfdot')
  res <- res@pcv
  rownames(res) <- colnames(exprData)
  return(res)
  
}
