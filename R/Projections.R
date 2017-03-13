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

performPCA <- function(data, N=0, variance_proportion=1.0) {
  #' Performs PCA on data
  #' 
  #' Parameters:
  #'  data: (Num_Features x Num_Samples) matrix
  #'    Matrix containing data to project into 2D
  #'  N: int
  #'    Number of Principle Components to reatin
  #'  variance_proportion: float
  #'    Retain top X principal components such taht a total of <variance_proportion> of the 
  #'    variance is retained  
  #'  
  #' Returns:
  #'  pca_data: (Num_Components x Num_Samples) matrix
  #'    Data transformed using PCA. Num_Components = Num_Samples
  
  datat = t(data)
  
  res <- prcomp(datat, center=TRUE, scale=TRUE)
  
  if(N == 0) {
    
    total_var <- as.matrix(cumsum(res$sdev^2 / sum(res$sdev^2)))
    last_i <- tail(which(total_var <= variance_proportion), n=1)
    N <- last_i
  }
  
  return (res$x[,1:N])
}
