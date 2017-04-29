## File containing scripts to convert FastProject data objects into JSON objects
## Called from methods-FastProjectOutput during the Output phase

require(jsonlite)


## Wrapper class for Expression data when generating json data
setMethod("initialize", signature(.Object="ServerExpression"), 
          function(.Object,data, sample_labels, gene_labels) {
            
            .Object@data <- data
            .Object@sample_labels <- sample_labels
            .Object@gene_labels <- gene_labels
            
            return(.Object)
          }
)

setMethod("initialize", signature(.Object="ServerSigProjMatrix"),
          function(.Object, data, proj_labels, sig_labels) {
            
            .Object@data <- data
            .Object@proj_labels <- proj_labels
            .Object@sig_labels <- sig_labels
            
            return(.Object)
          }
)

setMethod("initialize", signature(.Object="ServerPMatrix"),
          function(.Object, data, proj_labels, sig_labels) {
            
            .Object@data <- data
            .Object@proj_labels <- proj_labels
            .Object@sig_labels <- sig_labels
            
            return(.Object)
          })

signatureToJSON <- function(sig) {
  # Pass in a Signature object from a FastProjectOutput Object to be converted into a JSON object
  
  sig@sigDict <- as.list(sig@sigDict)
  
  json <- toJSON(sig, force=T, pretty=T, auto_unbox=T)
  return(json)
 
}

expressionToJSON <- function(expr, geneList=NULL) {
  #' Pass in an expression matrix from a FastProjectOutput Object
  #' Optionally, can specify a geneList which will select a subset of rows from the expression matrix
  
  if (!is.null(geneList)) {
    geneList = intersect(geneList, rownames(expr))
    expr <- expr[geneList,]
  }
  
  sExpr <- ServerExpression(unname(expr), colnames(expr), rownames(expr))
  
  ejson <- toJSON(sExpr, force=T, pretty=T, auto_unbox=T)

  return(ejson)
}

sigScoresToJSON <- function(ss) {
  #' Pass in a row of the SigMatrix from a FastProjectOutput Object
  #' By default will have colnames that can be set to list names while converting to a JSON object  

  s <- as.list(ss)
  names(s) <- rownames(as.matrix(ss))
  json <- toJSON(s, force=T, pretty=T, auto_unbox=T)
  
  return(json)
}

sigRanksToJSON <- function(ss) {
  #' Pass in a row of the SigMatrix from a FastProjectOuput object
  #' Compute ranks of sig scores from row and output JSON object mapping sample_label -> rank
  
  s <- as.list(rank(as.matrix(ss)))
  names(s) <- names(ss)
  json <- toJSON(s, force=T, pretty=T, auto_unbox=T)
  
  return(json)
  
}

coordinatesToJSON <- function(p) {
  #' Converts a projection into a JSON object mapping each sample to a projection coordinate
  
  coord <- apply(unname(p), 2, as.list)
  names(coord) <- colnames(p)
  
  json <- toJSON(coord, force=T, pretty=T, auto_unbox=T)
  
  return(json)
}

sigProjMatrixToJSON <- function(sigpm) {
  #' Converts a sigProjMatrix from an FastProjectOutput Object to a JSON object
  
  sSPM <- ServerSigProjMatrix(unname(sigpm), colnames(sigpm), rownames(sigpm))
  
  json <- toJSON(sSPM, force=T, pretty=T, auto_unbox=T)
  
  return(json)
}

sigProjMatrixPToJSON <- function(sigpmp) {
  #' Converts the -log10(pvalues) of the consistency scores into a JSON object
  
  sPM <- ServerPMatrix(unname(sigpmp), colnames(sigpmp), rownames(sigpmp))
  
  json <- toJSON(sPM, force=T, pretty=T, auto_unbox=T)
  
  return(json)

}

clusterToJSON <- function(cluster) {
  # Convert a cluster on some projection into a JSON object
  
  cl <- as.list(cluster)
  names(cl) <- colnames(cluster)
  
  json <- toJSON(cl, force=T, pretty=T, auto_unbox=T)
  
  return(json)
  
}

