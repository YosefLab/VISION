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
#' @param min_signature_genes Minimum number of genes to evaluate
#' @return SignatureScore object
naiveEvalSignature <- function(exprData, sig, weights, min_signature_genes) {

    # Select genes in signature that are in the expression matrix
    sig_names <- toupper(names(sig@sigDict))
    data_names <- toupper(rownames(exprData))
    sigVector <- sig@sigDict[sig_names %in% data_names]

    if (length(sigVector) == 0) {
    stop("No genes match signature.")
    }
    if (length(sigVector) < min_signature_genes) {
    stop("Too few genes match signature.")
    }

    sigGenes <- exprData[names(sigVector),]

    weights <- matrix(1, nrow=nrow(sigGenes), ncol=ncol(sigGenes))


    pdata <- sigGenes * sigVector * weights;

    sigScores <- colSums(pdata)

    sigScores <- sigScores / length(sigVector)

    sigObj <- SignatureScores(sigScores, sig@name, colnames(pdata),
                            isFactor=FALSE, isPrecomputed=FALSE, numGenes=length(sigVector))

    return(sigObj)
}

#' Evaluate signature with weights computed from FNR curve.
#'
#' @param exprData numeric Matrix Genes x Cells
#' @param sig Signature to be evalauting
#' @param weights Weight matrix computed through FNR curve
#' @param min_signature_genes Minimum number of genes to evaluate
#' @return SignatureScore object
weightedEvalSignature <- function(exprData, sig, weights, min_signature_genes) {

    # Select genes in signature that are in the expression matrix
    sig_names <- toupper(names(sig@sigDict))
    data_names <- toupper(rownames(exprData))
    sigVector <- (sig@sigDict)[sig_names %in% data_names]

    if (length(sigVector) == 0) {
    stop("No genes match signature.")
    }
    if (length(sigVector) < min_signature_genes) {
    stop("Too few genes match signature.")
    }

    sigGenes <- exprData[names(sigVector),]
    weights <- weights[names(sigVector),]

    pdata <- sigGenes * sigVector * weights

    sigScores <- colSums(pdata)
    denom <- colSums(abs(sigVector) * weights)
    denom[denom == 0] <- 1 # change 0/0 to 0/1

    sigScores <- sigScores / denom

    sigObj <- SignatureScores(sigScores, sig@name, colnames(pdata),
                            isFactor=FALSE, isPrecomputed=FALSE, numGenes=length(sigVector))

    return(sigObj)

}

## Define single signature evaluation for lapply method
singleSigEval <- function(s, sig_score_method, eData, weights, min_signature_genes) {
    # init to an empty SignatureScores object
        x <- NULL
        if (sig_score_method=="naive") {
            tryCatch({
                x <- naiveEvalSignature(eData, s, weights, min_signature_genes)
            }, error=function(e){})
        } else if (sig_score_method=="weighted_avg") {
            tryCatch({
                x <- weightedEvalSignature(eData, s, weights, min_signature_genes)
            }, error=function(e){})
        }
        return(x)
    }
