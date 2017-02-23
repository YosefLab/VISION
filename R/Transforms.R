## Functions to transform the data and calculate weights

createFalseNegativeMap <- function(data, housekeeping_genes) {
  #' Uses gene names in `housekeeping_genes` to create a mapping of false negatives.
  #' Creates a functional fit for each sample based on that samples HK genes
  #'
  #' Returns: fit_func: (function) Used to fit expression values to FN rate
  #'          params: (data.frame - Num_Params x Num_Samples) 
  #'                  Sample-specific parameters to use with fit_func
 
  # get subset of genes to be used, ie those included in the housekeeping genes set
  data_hk <- subset(data,  data[,1] %in% housekeeping_genes[[1]])[,-1]
  
  # convert data frame to numeric, and filter out genes with no variance
  data_hk <- sapply(data_hk, as.character)
  data_hk <- apply(data_hk, 2, as.numeric)
  data_hk <- filterGenesNovar(data_hk)
  
  # calculate the distributions for hk gene
  # Gamma is 1 for any non-zero data point
  # Mu_h is the row (per gene) average of non zero points
  gamma <- as.matrix(data_hk)
  gamma_i <- which(gamma > 0)
  gamma[gamma_i] <- 1
  mu_h <- as.matrix(apply(data_hk, 1, function(r) sum(r) / sum(r!=0)))
  
  
  # Fit a function mapping mu to gammas
  func <- function(xvals, x0, a, L=0, S=1) {
    return(L + (S/(1 + exp((xvals-x0)*a))))
  }
  
  
  efun <- function(x, y, args) {
    if (args[1] < 0) {
      args[1] = 0
    } else if (args[1] > Inf) {
      args[1] = Inf
    }
    
    if (args[2] < 0) {
      args[2] = 0
    } else if (args[2] > 2) {
      args[2] = 2
    }
    out <- func(x, args[1], args[2])
    return(sum((out-y)**2))
  }

  print(ncol(gamma))
  params <- matrix(0L, ncol=ncol(gamma), nrow=4)
  print(ncol(params))
  x <- c(mu_h)
  
  if(length(x) > 30) {
    q_indices <- round(length(x)/30 * seq(0, 29))
  } else {
    q_indices <- seq(0, 29)
  }

  q_indices <- c(q_indices, length(x))
  

  sort_i <- order(x)
  x_sorted <- x[sort_i]


  y <- 1-gamma
  y_sorted <- y[sort_i,]
  
  # Store the mean expression of genes per quantile
  x_quant <- rep(0, length(q_indices)-1);
  # Store the mean expression of genes in a sample per quantile
  y_quant <- matrix(0L, nrow=length(q_indices)-1, ncol=ncol(y))
  
  for(i in 1:(length(q_indices)-1)){
    start_i <- q_indices[i]+1;
    end_i <- q_indices[i+1];
    
    x_quant[i] <- mean(x_sorted[start_i:end_i]);
    y_quant[i,] = colMeans(y_sorted[start_i:end_i,])

  }
  
  bounds <- list(c(0, Inf), c(0, 2))
  initialGuesses <- list(c(3.5, 1), c(5.5, 1), c(1.5, .5), c(5.5, .5), c(3.5, 1.7))

  for (k in 1:(ncol(gamma) - 1)) {
    best_eval <- 1e99
    for (ig in initialGuesses) {
      
      res <- optim(par=c(ig), efun, x=x_quant, y=y_quant[,k])
      
      if (res$value < best_eval) {
        best_eval <- res$value
        print(best_eval)
        param <- res$par
        params[1,k] = param[1]
        params[2,k] = param[2]
        params[3,k] = 0
        params[4,k] = 1
        
      }
    }
  }
  
  return(c(func, params))
  
  
}