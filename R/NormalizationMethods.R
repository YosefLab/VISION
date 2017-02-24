#' Methods to normalize data
#' 
#' Different options to normalize the data before the signature
#' values are calculated.
#' 
#' Determined by the sig_norm_method value in the FastProject call
require(matrixStats)

noNormalization <- function(data) {
  #' Does nothing, return the original data
  
  message("No normalization applied")
  return(data)
}

colNormalization <- function(data) {
  #' Perform z-normalization on all columns
  
  message("Applying z-normalization on all columns...")

  # Normalize columns of data
  ndata <- scale(data)
  
  return(ndata)

}

rowNormalization <- function(data) {
  #' Perform z-normalization on all rows

  message("Applying z-normalization on all rows...")
  
  mu <- rowMeans(data)
  sigma <- rowSds(data)
  sigma[sigma == 0] <- 1.0
  
  ndata = ((data - mu) / sigma)
  
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