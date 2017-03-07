#' Functions for generating projections
#' 
#' This module handles the generation of lower-dimensional
#' projections from the higher-dimensional data objects.


generateProjections <- function(exprData, filterName="", inputProjections=NULL) {
  #' Projects data into 2 dimensions using a variety of linear and non-linear methods
  #' 
  #' Parameters:
  #'  exprData: (ExpressionData) expression matrix to project into 2D
  #'  filterName: (character) name of filter to apply to signatures, should match a filter that was
  #'              applied to the exprData before
  #'  inputProjections: (list of ProjectionData) collection of projection data types; names are of type
  #'              character, mapping to a projection             
  #' 
  #' Returns:
  #'  projections: (list of ProjectionData) dictionary mapping projection type to the actual projection
  #'  PC_data: (matrix) weighted PCA of original data type
  
  if (inputProjections == NULL) {
    inputProjections = c()
  }
  
  
  
  
  
}
