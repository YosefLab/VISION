#' Computes the biased SD on a vector, correcting the denominator to n rather than (n-1)
#' @importFrom stats sd
#' @param data Vector of numbers
#' @return Biased standard deviation of vector
biasedVectorSD <- function(data) {
    d <- length(data)
    std <- sd(data) * sqrt((d-1)/d)
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
#' @importFrom matrixStats colMeans2
#' @importFrom matrixStats colSds
#' @param data data matrix
#' @return Data matrix with same dimensions, with each column z-normalized.
colNormalization <- function(data) {

    message("Applying z-normalization on all columns...")

    mu <- colMeans2(data)
    sigma <- colSds(data)
    sigma[sigma == 0] <- 1.0
    ndata <- t(
               (t(data) - mu)/sigma
              )
    return(ndata)
}

#' Performs z-normalization on all rows
#' @importFrom matrixStats rowMeans2
#' @importFrom matrixStats rowSds
#' @param data data matrix
#' @return Data matrix with same dimensions, with each row z-normalized.
rowNormalization <- function(data) {

    message("Applying z-normalization on all rows...")
    mu <- rowMeans2(data)
    sigma <- rowSds(data)
    sigma[sigma == 0] <- 1.0
    ndata <- (data - mu) / sigma
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

    return(data)
}

#' Creaes a new version of the data that has ranks (column-wise) instead of values.
#' @importFrom matrixStats colRanks
#' @param data data matrix
#' @return Data matrix with same dimensions, with each value representing its column-wise rank.
colRankNormalization <- function(data) {

    message("Applying column rank normalization...")

    rdata = colRanks(data, ties.method="min",  preserveShape=TRUE)
    rownames(rdata) <- rownames(data)
    colnames(rdata) <- colnames(data)
    return(rdata)
}
