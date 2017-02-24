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
  return(data.fram(data))
}

colNormalization <- function(data) {
  #' Perform z-normalization on all columns
  
  message("Applying z-normalization on all columns...")
  
  ndata <- data[-1, -1]
  
  # Need to coerce elements of data frames for operations
  data <- apply(data, 2, as.character)
  ndata <- apply(ndata, 2, as.numeric)
  
  # Normalize columns of data
  ndata <- scale(ndata)
  #Coerce normalized data to insert back into orginal data file
  ndata <- apply(ndata,2, as.character)
 
  data[-1, -1] <- ndata
  return(data.frame(data))

}

rowNormalization <- function(data) {
  #' Perform z-normalization on all rows

  message("Applying z-normalization on all rows...")
  
  ndata <- data[-1,-1]
  
  data <- apply(data, 2, as.character)
  ndata <- apply(ndata, 2, as.numeric)
  
  mu <- rowMeans(ndata)
  sigma <- rowSds(ndata)
  sigma[sigma == 0] <- 1.0
  
  ndata = ((ndata - mu) / sigma)
  
  ndata <- apply(ndata, 2, as.character)
  data[-1, -1] <- ndata
  return(data.frame(data))
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
  
  ndata <- data[-1, -1]
  ndata <- apply(ndata, 2, as.numeric)
  data <- apply(data, 2, as.character)
  
  rdata <- matrix(0L, nrow=nrow(ndata), ncol=ncol(ndata))
  for (i in 1:ncol(rdata)) {
    rdata[,i] <- rank(ndata[,i], ties.method="min")
  }
  rdata <- apply(rdata, 2, as.character)
  data[-1,-1] <- rdata
  return(data.frame(data))
}