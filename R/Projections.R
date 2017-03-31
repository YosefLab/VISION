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


registerMethods <- function(lean=FALSE) {
  
  #projMethods <- c("ICA" = applyICA)
  projMethods <- c()
  if (!lean) {
    projMethods <- c(projMethods, "Spectral Embedding" = applySpectralEmbedding)
    projMethods <- c(projMethods, "MDS" = applyMDS)
  }
  
  #projMethods <- c(projMethods, "RBF Kernel PCA" = applyRBFPCA)
  #projMethods <- c(projMethods, "ISOMap" = applyISOMap)
  projMethods <- c(projMethods, "tSNE30" = applytSNE30)
  projMethods <- c(projMethods, "tSNE10" = applytSNE10)
  
  return(projMethods)
}

generateProjections <- function(expr, filterName="", inputProjections=c(), lean=FALSE) {
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
  
  print("PCA")
  pca_res <- applyPCA(exprData, N=30)
  proj <- Projection("PCA", pca_res)
  inputProjections <- c(inputProjections, proj)
  
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
    if (p@name == "PCA") {
      coordinates <- t(p@pData)
      coordinates <- coordinates[,1:2]
    } else {
      coordinates <- p@pData
    }
    
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

defineClusters <- function(projections) {
  #' Creates several different clusterings of the data in projections
  #' 
  #' Parameters:
  #'  projections: (list of Projection Objects) 
  #'    List of Projection objects, each of a different projection algorithm.
  #'  
  #' Returns:
  #'  clusters: (list of Cluster Objects)
  #'    list of cluster objects, each wrapping a projection with a groupings of different clusters

  # List of ClusterData Objects to be returned 
  outClusters <- c()
  
  for(proj in projections) {
    
    # List of Clusters for a given projection
    projClusters <- c()
    key <- proj@name
    
    for (k in 2:6) {
      clustName <- paste0("K-Means, k=", k)
      km <- kmeans(t(proj@pData), centers=k)
      clus <- Cluster(clustName, km$centers, t(as.matrix(km$cluster)))
      projClusters <- c(projClusters, clus)
    }
    
    cData <- ClusterData(key, projClusters)
    outClusters <- c(outClusters, cData)
  }
  
  return(outClusters)
  
  
  
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
  which(is.na(ndata), arr.ind=T)
  ndataT <- dimRedData(t(ndata))
  res <- embed(ndataT, "Isomap", knn=4, ndim=2)
  
  res <- res@data@data
  colnames(res) = colnames(exprData)
  return(t(res))
  
}

applyRBFPCA <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- kpca(ndataT, features=2, kernel='rbfdot')
  res <- res@pcv
  
  rownames(res) <- colnames(exprData)
  return(t(res))
  
}
