#' Calculates standard deviation with a denominator of n rather than (n-1)
#'
#' @param data data matrix to apply the biased standard deviation operation to
#' @param byRow if TRUE, apply by row. Else, by column. Default is TRUE.
#' @return Vector of standard deviations by row or column
biasedSD <- function(data, byRow=TRUE) {

  a <- 1
  d <- ncol(data)
  if (!byRow) {
    a = 2
    d <- nrow(data)
  }
  std <- apply(data, a, stats::sd) * sqrt(((d-1)/d))
  return(std)
}

#' Computes the biased SD on a vector, correcting the denominator to n rather than (n-1)
#'
#' @param data Vector of numbers
#' @return Biased standard deviation of vector
biasedVectorSD <- function(data) {
  d <- length(data)
  std <- stats::sd(data) * sqrt((d-1)/d)
  return(std)
}

#' Does nothing, just returns the original data
#'
#' @param data data matrix
#' @return The data matrix
noNormalization <- function(data) {
  return(data)
}

#' Performs z-normalization on all columns
#'
#' @param data data matrix
#' @return Data matrix with same dimensions, with each column z-normalized.
colNormalization <- function(data) {

  message("Applying z-normalization on all columns...")

  ndata <- t(data)
  mu <- rowMeans(ndata)
  sigma <- biasedSD(ndata)
  ndata <- (ndata - mu) / sigma
  return(t(ndata))
}

#' Performs z-normalization on all rows
#'
#' @param data data matrix
#' @return Data matrix with same dimensions, with each row z-normalized.
rowNormalization <- function(data) {

  message("Applying z-normalization on all rows...")
  mu <- rowMeans(data)
  sigma <- biasedSD(data)
  sigma[sigma == 0] <- 1.0
  ndata <- apply(data, 1, function(r) (r - mean(r)) / biasedVectorSD(r))
  ndata <- ( (data - mu) / sigma)
  return(ndata)
}

#' Performs z-normalization on all columns and rows
#'
#' @param data data matrix
#' @return Data matrix with same dimensions, with each row and column z-normalized.
rowAndColNormalization <- function(data) {
  message("Applying row and column normalization...")

  data = rowNormalization(data)
  data = colNormalization(data)

  return(data.frame(data))
}

#' Creaes a new version of the data that has ranks (column-wise) instead of values.
#'
#' @param data data matrix
#' @return Data matrix with same dimensions, with each value representing its column-wise rank.
colRankNormalization <- function(data) {

  message("Applying column rank normalization...")

  rdata <- matrix(0L, nrow=nrow(data), ncol=ncol(data))

  for (i in 1:ncol(rdata)) {
    rdata[,i] <- rank(data[,i], ties.method="min")
  }
  rownames(rdata) <- rownames(data)
  colnames(rdata) <- colnames(data)
  return(rdata)
}
