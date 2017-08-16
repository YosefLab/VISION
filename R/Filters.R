## Functions that are used to select genes

## Functions here use a variety of criteria to reduce the number
## of genes to a more manageable size - ideally extracting the
## genes that are more biologically informative.

#' Applies filters to the inputted expression data (may remove rows)
#'
#'  @param data ExpressionData object
#'  @param threshold minimum number of samples gene must be detected in to pass
#'  @param filterInput list of filters to compute
#'  @return List of the updated ExpressionData object, and the list of filters applied.
#'  @examples
#'  expr <- readExprToMatrix("data/expression.txt")
#'  edata <- ExpressionData(expr)
#'  filteredData <- applyFilters(eData, 10, c("threshold", "fano"))
applyFilters <- function(data, threshold, filterInput) {
  
  filterList <- c()
  expr <- getExprData(data)
  
  message("Applying filters...")
  
  for (filter in filterInput) {
    if (filter == "novar") {
      filterList <- c(filterList, "novar")
      expr <- filterGenesNovar(expr)
      data@noVarFilter <- expr
    } else if (filter == "threshold") {
      filterList <- c(filterList, "threshold")
      expr <- filterGenesThreshold(expr, threshold)
      data@thresholdFilter <- expr
    } else if (filter == "fano") {
      filterList <- c(filterList, "fano")
      expr <- filterGenesFano(expr)
      data@fanoFilter <- expr
    } else {
      stop("Filter not recognized")
    }
  }
  
  return(list(data, filterList))
  
}

#' Eliminate genes whose sample variance is equal to 0 (may remove rows); run when --nofilter option 
#' is selected
#' 
#' @param data expression matrix
#' @return filtered expression matrix
filterGenesNovar <- function(data) {

  message("Applying no variance filter...")
  

  return (subset(data, apply(data, 1, var) != 0))
  
}
#' Filter genes whose values sum to less than some threshold value (may remove rows)
#' 
#' @param data: (data.frame) expression matrix
#' @param threshold: (int) threshold value to filter by 
#' @return filtered expression matrix
filterGenesThreshold <- function(data, threshold) {

  message("Applying threshold filter...")
  
  return(subset(data, apply(data, 1, function(r) sum(r)) >= threshold))
  
  
}
# ' Applies the Fano filter to the input data (may remove rows)
#' 
#'	@param data NUM_GENES x NUM_SAMPLES expression matrix
#'	@param num_mad number of median absolute deviations 
#'  @return NUM_GENES_PASSED_FANO_FILTER x NUM_SAMPLES filtered expression matrix
filterGenesFano <- function(data, num_mad=2) {

  message("Applying fano filter...")
  sub_data <- data
  # if too many samples, subsample for fano filter
  if (ncol(data) > 50000) {
	sub_data <- data[,sample(ncol(data), 50000)]
  }
  
  mu <- apply(sub_data, 1, mean)
  sigma <- apply(sub_data, 1, sd)
  
  
  aa <- order(mu)
  mu_sort <- mu[aa]
  sigma_sort <- sigma[aa]
  
  
  N_QUANTS <- 30
  m <- floor(length(mu_sort) / N_QUANTS)
  
  gene_passes <- rep(0, nrow(sub_data)) == 1

  genePassList <- lapply(seq(0, N_QUANTS), function(i) {
		if (i == N_QUANTS - 1) {
			rr <- seq(i*m, length(mu_sort))
		} else {
			rr <- seq(i*m, (i+1)*m)
		}

		mu_quant <- mu_sort[rr]
		mu_quant[mu_quant==0] = 1
		sigma_quant <- sigma_sort[rr]
		fano_quant <- (sigma_quant ** 2) / (mu_quant)
		mad_quant <- median(abs(fano_quant - median(fano_quant)))
		gene_passes_quant <- (fano_quant > (median(fano_quant) + num_mad * mad_quant))
		gene_passes_quant_i <- which(gene_passes_quant != 0)
		gene_passes_i <- gene_passes_quant_i + (i*m)
		return(gene_passes_i)
	
	})

  gpi <- unlist(genePassList)
  gene_passes[gpi] <- T
  gene_passes <- gene_passes[order(aa)]
  return(data[which(gene_passes==T),])
  
}

