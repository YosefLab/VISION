## Functions to transform the data and calculate weights

probabilityOfExpression <- function(data) {
  
}

createFalseNegativeMap <- function(data, housekeeping_genes) {
  #' Uses gene names in `housekeeping_genes` to create a mapping of false negatives.
  #' Creates a functional fit for each sample based on that samples HK genes
  #'
  #' Returns: fit_func: (function) Used to fit expression values to FN rate
  #'          params: (data.frame - Num_Params x Num_Samples) 
  #'                  Sample-specific parameters to use with fit_func
  
  keep_indices <- c()
  
  for (hkgene in housekeeping_genes) {
    #for (gene in enumerate(data.row_labels)){
  
    #  if (lower(gene) == lower(hkgene)) {
    #    keep_indices <- c(keep_indices, i)
    #    continue;
    #  }
    
    #}
  }
  
}