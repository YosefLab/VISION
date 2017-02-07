## Functions that are used to select genes

## Functions here use a variety of criteria to reduce the number
## of genes to a more manageable size - ideally extracting the
## genes that are more biologically informative.


applyFilters <- function(data, threshold, nofilter, fano) {
  #'Applies filters to the inputted expression data (may remove rows)
  #'
  #'Parameters:
  #'  data: (data.frame) expression matrix
  #'  threshold: (int) minimum number of samples gene must be detected in to pass
  #'  nofilter: (logical) if true, only filter rows that have all identical values
  #'  fano: (logical) if true, apply the fano filter to the expression matrix
  #'Returns: (data.frame) filtered expression matrix
  
  
  if(nofilter) {
  
    data <- filterGenesNovar(data.frame(data))

  } else {

    if (!fano) {
      
      data <- filterGenesThreshold(data.frame(data), threshold)
      
    } else {
      
      data <- filterGenesFano(data.frame(data))
      
    }
    
  }
  
  col_labels <- data[1,]
  row_labels <- data[,1]
  
  return(list(data, row_labels, col_labels))
}

# Remove genes with 0 variance
filterGenesNovar <- function(data) {
  #' Eliminate genes whose sample variance is equal to 0 (may remove rows); run when --nofilter option 
  #' is selected
  #' 
  #' Parameters:
  #'  data: (data.frame) expression matrix
  #' Returns: (data.frame) filtered expression matrix
  
  d <- data.frame(data)
  return (subset(d, apply(d[,-1], 1, var) != 0))
  
}

filterGenesThreshold <- function(data, threshold) {
  #' Filter genes whose values sum to less than some threshold value (may remove rows)
  #' 
  #' Parameters:
  #'  data: (data.frame) expression matrix
  #'  threshold: (int) threshold value to filter by 
  #' Returns: (data.frame) filtered expression matrix

  d <- data.frame(data)
  return (subset(d, rowSums(d[,-1]) >= threshold))
  
}

filterGenesFano <- function(data, num_mad=2) {
  #' Applies the Fano filter to the input data (may remove rows)
  #' 
  #' Paramaters:
  #'  data: (data.frame) expression matrix
  #'  num_mad: (float) number of median absolute deviations 
  #' Returns: (data.frame) matrix with logical values, TRUE if gene passes fano filter; else FALSE
  
  d <- data.frame(data)
  
  mu <- apply(d[,-1], 1, mean)
  sigma <- apply(d[,-1], 1, sd)
  
  
  aa <- order(mu)
  mu_sort <- mu[aa]
  sigma_sort <- sigma[aa]
  
  
  N_QUANTS <- 30
  m <- floor(length(mu_sort) / N_QUANTS)
  
  gene_passes <- rep(0, nrow(d)) == 1
  
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
    gene_passes_quant_i = which(gene_passes_quant != 0);
    gene_passes_i <- gene_passes_quant_i + (i * m);
    gene_passes[gene_passes_i] <- TRUE

  }

  original_ii = order(aa)
  gene_passes <- gene_passes[original_ii]
  print(length(which(gene_passes == TRUE)))
  return(gene_passes)
}

