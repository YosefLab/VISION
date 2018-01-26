#' Different ways to evalutate the signature score
#'
#' EAch method should have the same signature so they can be swapped
#'
#' Right now, each takes in a wrapped data object and signature object
#'
#' Specified by FastProject argument (sig_score_method), default = naiveEvalSignature

#' Naive eval Siganture, just sums the columns * sign, equivalent to all weights = 1
#'
#' @param exprData numeric Matrix Genes x Cells
#' @param sig Signature to be evalauting
#' @param weights Weight matrix computed through FNR curve
#' @return SignatureScore object
naiveEvalSignature <- function(exprData, sig, weights) {

    # Select genes in signature that are in the expression matrix
    sig_names <- names(sig@sigDict)
    data_names <- rownames(exprData)
    sigVector <- sig@sigDict[sig_names %in% data_names]

    sigGenes <- exprData[names(sigVector),]

    weights <- matrix(1, nrow=nrow(sigGenes), ncol=ncol(sigGenes))


    pdata <- sigGenes * sigVector * weights;

    sigScores <- colSums(pdata)

    sigScores <- sigScores / length(sigVector)

    sigObj <- SignatureScores(sigScores, sig@name,
                            isFactor=FALSE, isPrecomputed=FALSE, numGenes=length(sigVector))

    return(sigObj)
}

#' Evaluate signature with weights computed from FNR curve.
#'
#' @param exprData numeric Matrix Genes x Cells
#' @param sig Signature to be evalauting
#' @param weights Weight matrix computed through FNR curve
#' @return SignatureScore object
weightedEvalSignature <- function(exprData, sig, weights) {

    # Select genes in signature that are in the expression matrix
    sig_names <- names(sig@sigDict)
    data_names <- rownames(exprData)
    sigVector <- (sig@sigDict)[sig_names %in% data_names]

    sigGenes <- exprData[names(sigVector),]
    weights <- weights[names(sigVector),]

    pdata <- sigGenes * sigVector * weights

    sigScores <- colSums(pdata)
    denom <- colSums(abs(sigVector) * weights)
    denom[denom == 0] <- 1 # change 0/0 to 0/1

    sigScores <- sigScores / denom

    sigObj <- SignatureScores(sigScores, sig@name,
                            isFactor=FALSE, isPrecomputed=FALSE, numGenes=length(sigVector))

    return(sigObj)

}

#' Evaluate signature scores efficiently in batches
#'
#' @param sigs list of Signature(s) to be evalauting
#' @param sig_score_method either "naive" or "weighted_avg"
#' @param eData numeric Matrix Genes x Cells
#' @param weights Weight matrix computed through FNR curve
#' @param min_signature_genes signatures are discarded which do not match at least
#' this many genes in the expression matrix
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @return List of SignatureScore objects
batchSigEval <- function(sigs, sig_score_method, eData, weights,
                         min_signature_genes) {

    workers <- BiocParallel::bpparam()$workers

    if (sig_score_method == "naive") {
        weights <- NULL
    }

    # Need to perform this multiply here so it doesn't occupy memory
    # in all of the sub-processes
    if( !is.null(weights) ) {
        expr_weights <- eData * weights
    } else {
        expr_weights <- eData
    }

    # Compute number of genes per signature that match and filter
    # the list of signatures

    numMatches <- vapply(sigs, function(sig){
        gene_count <- sum(names(sig@sigDict) %in% rownames(expr_weights))
        return(gene_count)
    }, 1.0)

    sigs <- sigs[numMatches >= min_signature_genes]

    # Partition signatures into batches
    # 1200 seems to be an ok batch size goal
    availableCores <- max(parallel::detectCores() - 1, 1)
    sigBatches <- batchify(sigs, 1200, n_workers = availableCores)

    allScoresBatches <- parallel::mclapply(sigBatches, function(sigBatch) {
        scores <- innerEvalSignatureBatch(expr_weights, sigBatch, weights)
        return(scores)
    }, mc.cores = min(availableCores, length(sigBatches)))

    # allScoresBatches is list of sig x cell matrices

    sigScoresObjBatches <- lapply(allScoresBatches, function(allScoresBatch) {

        sigScoresObj <- lapply(rownames(allScoresBatch), function(name) {
            sigObj <- SignatureScores(allScoresBatch[name, ], name,
                                    isFactor = FALSE, isPrecomputed = FALSE,
                                    numGenes = numMatches[name])
            return(sigObj)
        })

        return(sigScoresObj)
    })

    sigScores <- unlist(sigScoresObjBatches, recursive = FALSE, use.names = FALSE)

    names(sigScores) <- vapply(sigScores, function(s) s@name, "")

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
innerEvalSignatureBatch <- function(exprData, sigs, weights = NULL) {

    sigSparseMatrix <- sigsToSparseMatrix(sigs, exprData)

    sigScores <- sigSparseMatrix %*% exprData
    sigScores <- as.matrix(sigScores)

    if ( !is.null(weights) ) {
        denom <- abs(sigSparseMatrix) %*% weights # denom is N_sigs X N_cells
        denom <- as.matrix(denom)
        denom[denom == 0] <- 1
    } else  {
        denom <- rowSums(abs(sigSparseMatrix)) # denom is vector of length N_sigs
    }

    sigScores <- sigScores / denom

    return (sigScores)
}
