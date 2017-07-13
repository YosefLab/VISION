#' Different ways to evalutate the signature score
#' 
#' Each method should have the same signature so they can be swapped
#' 
#' Right now, each takes in a wrapped data object and signature object
#' 
#' Specified by FastProject argument (sig_score_method), default = weightedEvalSignature


naiveEvalSignature <- function(expr, sig, weights, min_signature_genes) {
  #' Naive eval signature, just sums the columns * sign
  #' Equivalent to all weights = 0
  #' 
  #' Paramters: exprData (ExpressionData) ExpressionData wrapper
  #'            sig (Signature) Sig to be evaluating
  #'            zeros (matrix) Locations of zeros in exprData
  #'            min_signature_genes (numeric) parameter set in FastProject, minimum number of genes to evaluate
  #'            
  #' Returns: sigObj (SignatureScore) Evaluation of signature
  
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

weightedEvalSignature <- function(expr, sig, weights, min_signature_genes) {
  #' Eval signature with weights stored in FastProject@weights
  #' 
  #' Paramters: exprData (ExpressionData) ExpressionData wrapper
  #'            sig (Signature) Sig to be evaluating
  #'            zeros (matrix) Locations of zeros in exprData
  #'            min_signature_genes (numeric) parameter set in FastProject, minimum number of genes to evaluate
  #'            
  #' Returns: sigObj (SignatureScore) Evaluation of signature
  
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
  weights <- weights[names(sigVector),]

  ## UNCOMMENT OUT THIS LINE TO ENABLE C++ SIG SCORE CALCULATION
  #pdata <- calcSigScore(sigGenes, sigVector, weights)
  pdata <- sigGenes * sigVector * weights

  sigScores <- colSums(pdata)
  sigScores <- sigScores / colSums(abs(sigVector) * weights)
  
  sigObj <- SignatureScores(sigScores, sig@name, colnames(weights), 
                            isFactor=FALSE, isPrecomputed=FALSE, numGenes=length(sigVector))
  
  return(sigObj)
  
}
