## Functions that are used to select genes

## Functions here use a variety of criteria to reduce the number
## of genes to a more manageable size - ideally extracting the
## genes that are more biologically informative.

applyFilters <- function(data, threshold, nofilter, lean) {
  #'Applies filters to the inputted expression data (may remove rows)
  #'
  #'Parameters:
  #'  data: (ExpressionData) expression matrix
  #'  threshold: (int) minimum number of samples gene must be detected in to pass
  #'  nofilter: (logical) if true, only filter rows that have all identical values
  #'  lean: (logical) if true, skip extra filtering methods (Fano)
  #'Returns: (data.frame) filtered expression matrix
  
  expr <- getExprData(data)
  
  message("Applying filters...")
  
  if(nofilter) {
    expr <- filterGenesNovar(expr)
    data@noVarFilter <- expr
  } else {

    if (lean) {
      expr <- filterGenesThreshold(expr, threshold)
      data@thresholdFilter <- expr
    } else {
      expr <- filterGenesThreshold(expr, threshold)
      data@thresholdFilter <- expr
      
      expr <- filterGenesFano(expr)
      data@fanoFilter <- expr
    }
    
  }

  return(expr)
  
}


filterGenesNovar <- function(data) {
  #' Eliminate genes whose sample variance is equal to 0 (may remove rows); run when --nofilter option 
  #' is selected
  #' 
  #' Parameters:
  #'  data: (data.frame) expression matrix
  #' Returns: (data.frame) filtered expression matrix
  message("Applying no variance filter...")
  
  #d <- data.frame(data)
  return (subset(data, apply(data, 1, var) != 0))
  
}

filterGenesThreshold <- function(data, threshold) {
  #' Filter genes whose values sum to less than some threshold value (may remove rows)
  #' 
  #' Parameters:
  #'  data: (data.frame) expression matrix
  #'  threshold: (int) threshold value to filter by 
  #' Returns: (data.frame) filtered expression matrix
  #' 
  #' 
  message("Applying threshold filter...")
  
  return(subset(data, apply(data, 1, function(r) sum(r)) >= threshold))
  
  
}

filterGenesFano <- function(data, num_mad=2) {
  #' Applies the Fano filter to the input data (may remove rows)
  #' 
  #' Paramaters:
  #'  data: (data.frame) expression matrix
  #'  num_mad: (float) number of median absolute deviations 
  #' Returns: (data.frame) matrix with logical values, TRUE if gene passes fano filter; else FALSE
  message("Applying fano filter...")
  
  mu <- apply(data, 1, mean)
  sigma <- apply(data, 1, sd)
  
  
  aa <- order(mu)
  mu_sort <- mu[aa]
  sigma_sort <- sigma[aa]
  
  
  N_QUANTS <- 30
  m <- floor(length(mu_sort) / N_QUANTS)
  
  gene_passes <- rep(0, nrow(data)) == 1
  
  for (i in 0:N_QUANTS) {
    if (i == N_QUANTS - 1) {
      rr <- seq(from=i*m, to=length(mu_sort))
    } else {
      rr <- seq(i*m, (i+1) * m)
    }
    
    mu_quant <- mu_sort[rr]
    mu_quant[mu_quant==0] = 1
    sigma_quant <- sigma_sort[rr]
    fano_quant <- (sigma_quant ** 2) / (mu_quant)
    mad_quant <- median(abs(fano_quant - median(fano_quant)))
    gene_passes_quant <- fano_quant > (median(fano_quant) + num_mad * mad_quant)
    gene_passes_quant_i = which(gene_passes_quant != 0)
    gene_passes_i <- gene_passes_quant_i + (i * m);
    gene_passes[gene_passes_i] <- TRUE

  }

  original_ii = order(aa)
  gene_passes <- gene_passes[original_ii]
  print(length(which(gene_passes == TRUE)))
  keep_ii <- which(gene_passes == TRUE)
  return(data[keep_ii,])
}

