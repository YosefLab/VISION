#' Evaluate signature scores efficiently in batches
#'
#' This version uses the NormData object to operate
#' without inflating sparse matrices
#'
#' @param sigs list of Signature(s) to be evalauting
#' @param normData NormData object
#' @importFrom pbmcapply pbmclapply
#' @return matrix of signature scores, cells X signatures
batchSigEvalNorm <- function(sigs, normData) {

    # Partition signatures into batches
    # 1200 seems to be an ok batch size goal
    n_workers <- getOption("mc.cores")
    n_workers <- if (is.null(n_workers)) 2 else n_workers
    sigBatches <- batchify(sigs, 1200, n_workers = n_workers)

    allScoresBatches <- pbmclapply(sigBatches, function(sigBatch) {
        scores <- innerEvalSignatureBatchNorm(normData, sigs)
        return(scores)
    })

    # allScoresBatches is list of sig x cell matrices
    sigScores <- t(do.call(rbind, allScoresBatches))

    return(sigScores)
}

#' Utility method to load signatures into a sparse matrix
#'
#' @importFrom Matrix sparseMatrix
#' @param sigs List of Signature to be evalauting
#' @param expression numeric Matrix Genes x Cells
#' @return CsparseMatrix containing signature values
sigsToSparseMatrix <- function(sigs, expression) {

    sigMatches <- lapply(seq(length(sigs)), function(i) {
        sig <- sigs[[i]]
        indices <- match(names(sig@sigDict), rownames(expression))
        values <- unname(sig@sigDict)
        valid <- !is.na(indices)
        indices <- indices[valid]
        values <- values[valid]
        ii <- rep(i, length(indices))
        return(list(indices, ii, values))
    })

    j <- unlist(lapply(sigMatches, function(x) x[[1]]))
    i <- unlist(lapply(sigMatches, function(x) x[[2]]))
    x <- unlist(lapply(sigMatches, function(x) x[[3]]))

    dn <- list( vapply(sigs, function(x) x@name, ""),
                rownames(expression))

    sigSparseMatrix <- sparseMatrix(i = i, j = j, x = x,
                                    dims = c(length(sigs),
                                    nrow(expression)),
                                    dimnames = dn)

    return(sigSparseMatrix)
}


#' Used in inner loop of batchSigEvalNorm
#'
#' Computes signature scores without inflating the genes/cells matrix
#'
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#'
#' @param normData NormData row/column normalization factors
#' @param sigs List of Signature to be evalauting
#' @return matrix containing signature values (sigs x cells)
innerEvalSignatureBatchNorm <- function(normData, sigs) {

    sigSparseMatrix <- sigsToSparseMatrix(sigs, normData@data)

    NCells <- ncol(normData@data)
    NGenes <- nrow(normData@data)
    Rs <- Diagonal(x = normData@rowScaleFactors)
    Cs <- Diagonal(x = normData@colScaleFactors)
    Rog <- Matrix(normData@rowOffsets, ncol = 1)
    Roc <- Matrix(1, nrow = 1, ncol = NCells)

    SRs <- (sigSparseMatrix %*% Rs)
    SRsE <- SRs %*% normData@data
    SRsRo <- (SRs %*% Rog) %*% Roc

    # Note: this requires sparse=TRUE so OMP/MKL won't use many threads
    # for the next multiply (e.g., avoid multiplying dense x dense).  This
    # is important because this runs inside a parallel loop already
    Cog <- Matrix(1, ncol = 1, nrow = NGenes, sparse = TRUE)
    Coc <- Matrix(normData@colOffsets, nrow = 1, sparse = TRUE)

    SCo <- (sigSparseMatrix %*% Cog) %*% Coc

    C <- (SRsE + SRsRo + SCo) %*% Cs

    sigScores <- as.matrix(C)

    denom <- rowSums(abs(sigSparseMatrix)) # denom is vector of length N_sigs

    sigScores <- sigScores / denom

    return (sigScores)
}
