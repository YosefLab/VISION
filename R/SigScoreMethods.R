#' Different ways to evalutate the signature score
#'
#' EAch method should have the same signature so they can be swapped
#'
#' Right now, each takes in a wrapped data object and signature object
#'
#' Specified by FastProject argument (sig_score_method), default = naiveEvalSignature

#' Naive eval Siganture, just sums the columns * sign, equivalent to all weights = 1
#'
#' @param expr ExpressionData object
#' @param sig Signature to be evalauting
#' @param weights Weight matrix computed through FNR curve
#' @param min_signature_genes Minimum number of genes to evaluate
#' @return SignatureScore object
naiveEvalSignature <- function(expr, sig, weights, min_signature_genes) {

  exprData <- getExprData(expr)

  # Select genes in signature that are in the expression matrix
  keep_ii <- which(names(sig@sigDict) %in% rownames(exprData))
  sigVector <- (sig@sigDict)[keep_ii]

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
#' @param expr ExpressionData object
#' @param sig Signature to be evalauting
#' @param weights Weight matrix computed through FNR curve
#' @param min_signature_genes Minimum number of genes to evaluate
#' @return SignatureScore object
weightedEvalSignature <- function(expr, sig, weights, min_signature_genes) {

  exprData <- getExprData(expr)

  # Select genes in signature that are in the expression matrix
  keep_ii <- names(sig@sigDict) %in% rownames(exprData)
  sigVector <- (sig@sigDict)[keep_ii]

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
  sigScores <- sigScores / colSums(abs(sigVector) * weights)

  sigObj <- SignatureScores(sigScores, sig@name, colnames(pdata),
                            isFactor=FALSE, isPrecomputed=FALSE, numGenes=length(sigVector))

  return(sigObj)

}
