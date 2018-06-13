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

    data = rowNormalization(data)
    data = colNormalization(data)

    return(data)
}

#' Creaes a new version of the data that has ranks (column-wise) instead of values.
#' @importFrom matrixStats colRanks
#' @param data data matrix
#' @return Data matrix with same dimensions, with each value representing its column-wise rank.
colRankNormalization <- function(data) {

    rdata = colRanks(data, ties.method="min",  preserveShape=TRUE)
    rownames(rdata) <- rownames(data)
    colnames(rdata) <- colnames(data)
    return(rdata)
}

#' Calculates the specified normalized data matrix
#' @param data numeric matrix
#' @param func normalization method to apply
#' @return Normalized data matrix according to function specified.
getNormalizedCopy <- function(data, func) {

    data <- matLog2(data)

    if (func == "none") {
        return(noNormalization(data))
    } else if (func == "znorm_columns") {
        return(colNormalization(data))
    } else if (func == "znorm_rows") {
        return(rowNormalization(data))
    } else if (func == "znorm_rows_then_columns") {
        return(rowAndColNormalization(data))
    } else if (func == "rank_norm_columns") {
        return(colRankNormalization(data))
    }
    stop("Normalization method not recognized.")
}
