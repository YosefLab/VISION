#' Different ways to evalutate the signature score
#'
#' EAch method should have the same signature so they can be swapped
#'
#' Right now, each takes in a wrapped data object and signature object
#'
#' Specified by FastProject argument (sig_score_method), default = naiveEvalSignature

#' Evaluate signature scores efficiently in batches
#'
#' @param sigs list of Signature(s) to be evalauting
#' @param sig_score_method either "naive" or "weighted_avg"
#' @param eData numeric Matrix Genes x Cells
#' @param weights Weight matrix computed through FNR curve
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @return matrix of signature scores, cells X signatures
batchSigEval <- function(sigs, sig_score_method, eData, weights) {

    workers <- BiocParallel::bpparam()$workers

    if (sig_score_method == "naive") {
        weights <- matrix(NA, 1, 1)
    }

    # Need to perform this multiply here so it doesn't occupy memory
    # in all of the sub-processes
    if (!all(dim(weights) == c(1, 1))) {
        expr_weights <- eData * weights
    } else {
        expr_weights <- eData
    }

    # Partition signatures into batches
    # 1200 seems to be an ok batch size goal
    availableCores <- min(max(parallel::detectCores() - 1, 1), 10)
    sigBatches <- batchify(sigs, 1200, n_workers = availableCores)

    allScoresBatches <- parallel::mclapply(sigBatches, function(sigBatch) {
        scores <- innerEvalSignatureBatch(expr_weights, sigBatch, weights)
        return(scores)
    }, mc.cores = min(availableCores, length(sigBatches)))

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


#' Used in inner loop of batchSigEval
#'
#' @param exprData numeric Matrix Genes x Cells
#' @param sigs List of Signature to be evalauting
#' @param weights numeric Matrix Genes x Cells
#' @return matrix containing signature values (sigs x cells)
innerEvalSignatureBatch <- function(exprData, sigs, weights = matrix(NA, 1, 1)) {

    sigSparseMatrix <- sigsToSparseMatrix(sigs, exprData)

    sigScores <- sigSparseMatrix %*% exprData
    sigScores <- as.matrix(sigScores)

    if (!all(dim(weights) == c(1, 1))) {
        denom <- abs(sigSparseMatrix) %*% weights # denom is N_sigs X N_cells
        denom <- as.matrix(denom)
        denom[denom == 0] <- 1
    } else  {
        denom <- rowSums(abs(sigSparseMatrix)) # denom is vector of length N_sigs
    }

    sigScores <- sigScores / denom

    return (sigScores)
}
