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

#' Calculates the specified normalized data matrix
#'
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#'
#' @param data numeric matrix
#' @param func normalization method to apply
#' @return NormData normalized data object that maintains sparse structure
getNormalizedCopySparse <- function(data, func) {

    data <- matLog2(data)

    if (!(func %in% c("none", "znorm_columns", "znorm_rows",
        "znorm_rows_then_columns", "rank_norm_columns"))){
        stop("Normalization method not recognized.")
    }

    rowOffsets <- NULL
    colOffsets <- NULL
    rowScaleFactors <- NULL
    colScaleFactors <- NULL

    if (func == "znorm_columns" || func == "znorm_rows_then_columns") {
        colOffsets <- colMeans(data)
        colScaleFactors <- colVarsSp(data) ** 0.5
    }

    if (func == "znorm_rows" || func == "znorm_rows_then_columns") {
        rowOffsets <- rowMeans(data)
        rowScaleFactors <- rowVarsSp(data) ** 0.5
    }

    nd <- NormData(data, rowOffsets = rowOffsets, colOffsets = colOffsets,
        rowScaleFactors = rowScaleFactors, colScaleFactors = colScaleFactors)

    return(nd)
}

#' Initialize a new NormData object
#'
#' @param data expression data matrix
#' @param rowOffsets offsets to be subtracted from each row
#' @param rowScaleFactors factors to scale each row
#' @param colOffsets offsets to be subtracted from each column
#' @param colScaleFactors factors to scale each column
#' @return NormData object
NormData <- function(data, rowOffsets = NULL, rowScaleFactors = NULL,
    colOffsets = NULL, colScaleFactors = NULL) {

    if (is.null(rowOffsets)){
        rowOffsets <- numeric(nrow(data))
    }

    if (is.null(rowScaleFactors)){
        rowScaleFactors <- numeric(nrow(data)) + 1
    }

    if (is.null(colOffsets)){
        colOffsets <- numeric(ncol(data))
    }

    if (is.null(colScaleFactors)){
        colScaleFactors <- numeric(ncol(data)) + 1
    }

    .Object <- new("NormData",
        data = data,
        rowOffsets = rowOffsets,
        colOffsets = colOffsets,
        rowScaleFactors = rowScaleFactors,
        colScaleFactors = colScaleFactors
    )

    return(.Object)
}
