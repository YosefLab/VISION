#' Methods to normalize data
#' 
#' Different options to normalize the data before the signature
#' values are calculated.
#' 
#' Determined by the sig_norm_method value in the FastProject call
require(matrixStats)

biasedSD <- function(data, byRow=TRUE) {
  a <- 1
  d <- ncol(data)
  if (!byRow) {
    a = 2
    d <- nrow(data)
  }
  std <- apply(data, a, sd) * sqrt(((d-1)/d))
  return(std)
}

biasedVectorSD <- function(data) {
  d <- length(data)
  std <- sd(data) * sqrt((d-1)/d)
  return(std)
}


noNormalization <- function(data) {
  #' Does nothing, return the original data
  
  message("No normalization applied")
  return(data)
}

colNormalization <- function(data) {
  #' Perform z-normalization on all columns
  
  message("Applying z-normalization on all columns...")

  #ndata <- apply(data, 2, function(c) (c - mean(c)) / biasedVectorSD(c))
  #return(ndata)
  ndata <- t(data)
  mu <- apply(ndata, 1, mean)
  sigma <- biasedSD(ndata)
  ndata <- (ndata - mu) / sigma
  return(t(ndata))
}

rowNormalization <- function(data) {
  #' Perform z-normalization on all rows

  message("Applying z-normalization on all rows...")
  mu <- rowMeans(data)
  sigma <- biasedSD(data)
  sigma[sigma == 0] <- 1.0
  ndata <- apply(data, 1, function(r) (r - mean(r)) / biasedVectorSD(r))
  #return(ndata)
  ndata <- ( (data - mu) / sigma)
  return(ndata)
}

rowAndColNormalization <- function(data) {
  #' Normalize rows, then columns
  
  message("Applying row and column normalization...")
  
  data = rowNormalization(data)
  data = colNormalization(data)
  
  return(data.frame(data))
}

colRankNormalization <- function(data) {
  #' Create new version of data that has ranks (column-wise) instead of values

  message("Applying column rank normalization...")

  
  rdata <- matrix(0L, nrow=nrow(data), ncol=ncol(data))
  
  for (i in 1:ncol(rdata)) {
    rdata[,i] <- rank(data[,i], ties.method="min")
  }
  rownames(rdata) <- rownames(data)
  colnames(rdata) <- colnames(data)
  return(rdata)
}