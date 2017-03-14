#' Functions for generating projections
#' 
#' This module handles the generation of lower-dimensional
#' projections from the higher-dimensional data objects.

require("fastICA")
require('tsne')
require('dimRed')
require("RANN")
require('igraph')
require('kernlab')


registerMethods <- function(lean=FALSE) {
  
  projMethods <- c("ICA" = applyICA)
  
  if (!lean) {
    projMethods <- c(projMethods, "Spectral Embedding" = applySpectralEmbedding)
    projMethods <- c(projMethods, "MDS" = applyMDS)
  }
  
  projMethods <- c(projMethods, "RBF Kernel PCA" = applyRBFPCA)
  projMethods <- c(projMethods, "ISOMap" = applyISOMap)
  projMethods <- c(projMethods, "tSNE30" = applytSNE30)
  projMethods <- c(projMethods, "tSNE10" = applytSNE10)
  projMethods <- c(projMethods, "PCA" = applyPCA)
  
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
  
  for (method in names(methodList)){
    
    # Check if method is PCA because function call isn't general
    if (method == "PCA") {
      pca_res <- applyPCA(exprData, N=3)
      proj <- Projection("PCA", pca_res)
      inputProjections <- c(inputProjections, proj)
    }
    else {
      res <- methodList[[method]](exprData)
      proj <- Projection(method, res)
      inputProjections <- c(inputProjections, proj)
    }
  }  
  
  
  return(list(inputProjections, rownames(exprData)))
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
      clustName = paste0("K-Means, k=", k)
      km <- kmeans(t(proj@pData), centers=k)
      clus <- Cluster(clustName, km$centers, t(as.matrix(km$cluster)))
      projClusters <- c(projClusters, clus)
    }
    
    cData <- ClusterData(key, projClusters)
    outClusters <- c(outClusters, cData)
  }
  
  return(outClusters)
  
  
  
}

applyPCA <- function(data, N=0, variance_proportion=1.0) {
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
  
  datat = t(data)
  
  res <- prcomp(datat, center=TRUE, scale=TRUE)
  
  if(N == 0) {
    
    total_var <- as.matrix(cumsum(res$sdev^2 / sum(res$sdev^2)))
    last_i <- tail(which(total_var <= variance_proportion), n=1)
    N <- last_i
  }
  
  return (t(res$x[,1:N]))
}

applyICA <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- fastICA(ndataT, n.comp=2, maxit=100, tol=.0001, alg.typ="parallel", fun="logcosh", alpha=1,
                 method = "R", row.norm=FALSE, verbose=TRUE)
  
  return(t(res$S))
}

applySpectralEmbedding <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- specc(ndataT, centers=2)
  res <- specc(ndata, centers=2)
  res <- as.matrix(res@centers)
  
  colnames(res) = colnames(exprData)
  
  return(res)
  
}

applyMDS <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  distN <- dist(ndataT)
  
  res <- cmdscale(distN, k=2)
  return(t(res))
}

applytSNE10 <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- tsne(ndataT, k=2, perplexity=10.0, max_iter=200)
  
  rownames(res) = colnames(exprData)
  return(t(res))
  
}

applytSNE30 <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- t(ndata)
  res <- tsne(ndataT, k=2, perplexity=30.0, max_iter=200)
  
  rownames(res) = colnames(exprData)
  
  return(t(res))
}

applyISOMap <- function(exprData, projWeights=NULL) {
  set.seed(RANDOM_SEED)
  
  ndata <- colNormalization(exprData)
  ndataT <- dimRedData(t(ndata))
  res <- embed(ndataT, "Isomap", knn=4, ndim=2)
  
  res <- res@data@data
  rownames(res) = colnames(exprData)
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
