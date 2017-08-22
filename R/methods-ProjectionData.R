#' Wrapper for storing all relevant information for a given projection.
#' 
#' Stores a list of Projection objects, filter name, and a logical value indicating whether or not 
#' PCA was performed. Also stores clusters, a signature to projection matrix, and relevant gene names
#' and signature / projection keys.

#' Initializes a ProjectionData object for neatly storing all relevant data to a given projection section
#' 
#' @param filter Name of the filter applied to data before computing projections
#' @param projections List of Projection objects to be stored
#' @param genes Genes used to compute the projections
#' @param keys Sample names of expression data
#' @param sigProjMatrix Matrix storing the median consistency score per projection, signature pair
#' @param pMatrix Matrix storing the p values for each projection, signature pair
#' @param PPT List of Simple PPT parameters
#' @return ProjectionData object
setMethod("initialize", signature(.Object="ProjectionData"),
          function(.Object, filter = "", projections=NULL, genes, keys, sigProjMatrix, 
					pMatrix, PPT, fullPCA, mutualInformation, loadings) {
            
            .Object@filter <- filter
            .Object@projections <- projections
            .Object@genes <- genes
            .Object@keys <- keys
            .Object@sigProjMatrix <- sigProjMatrix
            .Object@pMatrix <- pMatrix
            .Object@PPT <- PPT
            .Object@fullPCA <- fullPCA
            .Object@mutualInformation <- mutualInformation
            .Object@loadings <- loadings
            
            return(.Object)
            
            
          }
)

