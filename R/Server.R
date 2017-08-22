## File containing scripts to convert FastProject data objects into JSON objects
## Called from methods-FastProjectOutput during the Output phase
	
require(jsonlite)

#' Wrapper class for ExpressionData object for JSON.
#' 
#' @param data Expression data 
#' @param sample_labels Labels of samples in expression data
#' @param gene_labels Lables of genes in expression data
#' @return ServerExpression object
setMethod("initialize", signature(.Object="ServerExpression"), 
          function(.Object,data, sample_labels, gene_labels) {
            
            .Object@data <- data
            .Object@sample_labels <- sample_labels
            .Object@gene_labels <- gene_labels
            
            return(.Object)
          }
)

#' Wrapper class for Signature Projection Matrix
#' 
#' @param data Signatrue Projection matrix as obtained from sigsVsProjections
#' @param proj_labels Projection names
#' @param sig_labels Names of signatures
#' @return ServerSigProjMatrix object
setMethod("initialize", signature(.Object="ServerSigProjMatrix"),
          function(.Object, data, proj_labels, sig_labels) {
            
            .Object@data <- data
            .Object@proj_labels <- proj_labels
            .Object@sig_labels <- sig_labels
            
            return(.Object)
          }
)

#' Wrapper class for the P value matrix calculated during sigsVsProjections
#' 
#' @param data P values for each signature, projection pair in the form of a matrix
#' @param proj_labels Projection names
#' @param sig_labels Names of signatures
#' @return ServerPMatrix object
setMethod("initialize", signature(.Object="ServerPMatrix"),
          function(.Object, data, proj_labels, sig_labels) {
            
            .Object@data <- data
            .Object@proj_labels <- proj_labels
            .Object@sig_labels <- sig_labels
            
            return(.Object)
          })


setMethod("initialize", signature(.Object="ServerMI"),
		function(.Object, data, proj_labels, sig_labels) {
		
		.Object@data <- data
		.Object@proj_labels <- proj_labels
		.Object@sig_labels <- sig_labels

		return(.Object)

		})


#' Converts Signature object to JSON
#' 
#' @param sig Signature object
#' @return JSON formatted Signature object.
signatureToJSON <- function(sig) {

  # Pass in a Signature object from a FastProjectOutput Object to be converted into a JSON object
  sig@sigDict <- as.list(sig@sigDict)
  
  json <- toJSON(sig, force=T, pretty=T, auto_unbox=T)
  return(json)
 
}

#' Convertes expression matrix to JSON
#' 
#' @param expr Expression Matrix
#' @param geneList optional list of genes to subset from expr
#' @return (Potentially) subsetted expression matrix
expressionToJSON <- function(expr, geneList=NULL) {

  if (!is.null(geneList)) {
    geneList = intersect(geneList, rownames(expr))
    expr <- expr[geneList,]
  }
  
  sExpr <- ServerExpression(unname(expr), colnames(expr), rownames(expr))
  
  ejson <- toJSON(sExpr, force=T, pretty=T, auto_unbox=T)

  return(ejson)
}

#' Converts row of sigantures score matrix to JSON
#' 
#' @param ss List of signature scores
#' @return Signature scores list to JSON, with names of each entry that of the list names
sigScoresToJSON <- function(ss) {

  s <- as.list(ss)
  names(s) <- rownames(as.matrix(ss))
  json <- toJSON(s, force=T, pretty=T, auto_unbox=T)
  
  return(json)
}

#' Converts list of signature ranks to JSON
#' 
#' @param ss List of rank values
#' @return Signature ranks as JSON, with names of each entry that of list names
sigRanksToJSON <- function(ss) {

  s <- as.list(rank(as.matrix(ss)))
  names(s) <- names(ss)
  json <- toJSON(s, force=T, pretty=T, auto_unbox=T)
  
  return(json)
  
}

#' Converts a projection into a JSON object mapping each sample to a projection coordinate.
#' 
#' @param p: Projection coordinate data (NUM_SAMPLES x NUM_COMPONENTS)
#' @return JSON object mapping each sample to a projection coordinate. 
coordinatesToJSON <- function(p) {
  
  coord <- apply(unname(p), 2, as.list)
  names(coord) <- colnames(p)
  
  json <- toJSON(coord, force=T, pretty=T, auto_unbox=T)
  
  return(json)
}

#' Converts a sigProjMatrix from a FastProjectOutput Object to a JSON object
#' 
#' @param sigpm SigProjMatrix
#' @param sigs Signatures to subset form sigpm
#' @return Subsetted sigProjMatirx converted to JSON
sigProjMatrixToJSON <- function(sigpm, sigs) {
  
  sigpm <- sigpm[sigs,, drop=F]
  sSPM <- ServerSigProjMatrix(unname(sigpm), colnames(sigpm), sigs)
  
  json <- toJSON(sSPM, force=T, pretty=T, auto_unbox=T)
  
  return(json)
}

mutualInfoToJSON <- function(mI, sigs) {
	
	mI <- mI[sigs,,drop=F]
	cn <- c()
	for (i in 1:ncol(mI)) { cn <- c(cn, paste("PC", i)) } 
	sMI <- ServerMI(unname(mI), cn, sigs)

	json <- toJSON(sMI, force=T, pretty=T, auto_unbox=T)

	return(json)

}

#' Converts the -log10(pvalues) of the consistency scores into a JSON object
#' 
#' @param sigpmp SigProjMatrix p values
#' @param sigs Signatrues to subset from sigpmp
#' @return Subsetted sigProjMatrix_P converted to JSON
sigProjMatrixPToJSON <- function(sigpmp, sigs) {

  sigpmp <- as.matrix(sigpmp[sigs,, drop=F])
  sPM <- ServerPMatrix(unname(sigpmp), colnames(sigpmp), sigs)
  
  json <- toJSON(sPM, force=T, pretty=T, auto_unbox=T)
  
  return(json)

}

#' Convert a Cluster object to JSON
#' 
#' @param cluster Cluster object
#' @return Cluster object converted to JSON
clusterToJSON <- function(cluster) {

  out <- list()
  out[['method']] <- cluster@method
  out[['param']] <- cluster@param
  out[['centers']] <- cluster@centers
  out[['data']] <- as.list(cluster@data[1,])
  json <- toJSON(out, force=T, pretty=T, auto_unbox=TRUE)

  return(json)
}

newAnalysis <- function(nfp) {
	fpo1 <- Analyze(nfp)
	saveFPOutAndViewResults(fpo1)

}
