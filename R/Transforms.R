## Functions to transform the data and calculate weights

createFalseNegativeMap <- function(data, housekeeping_genes) {
  #' Uses gene names in `housekeeping_genes` to create a mapping of false negatives.
  #' Creates a functional fit for each sample based on that samples HK genes
  #'
  #' Returns: fit_func: (function) Used to fit expression values to FN rate
  #'          params: (data.frame - Num_Params x Num_Samples) 
  #'                  Sample-specific parameters to use with fit_func
 
  # get subset of genes to be used
  
  data_hk <- subset(data,  data[,1] %in% housekeeping_genes[[1]])[,-1]
  
  data_hk <- sapply(data_hk, as.character)
  data_hk <- apply(data_hk, 2, as.numeric)
  data_hk <- filterGenesNovar(data_hk)
  
  gamma <- as.matrix(data_hk)
  gamma_i <- which(gamma > 0)
  gamma[gamma_i] <- 1
  mu_h <- as.matrix(apply(data_hk, 1, function(r) sum(r) / sum(r!=0)))

  
  
  func <- function(xvals, x0, a, L=0, S=1) {
    return(L + (S/(1 + exp((xvals-x0)*a))))
  }
  
  efun <- function(x, y, args) {
    out <- func(x, args[0], args[1])
    return(sum((out-y)**2))
  }

  
  params <- matrix( rep(0, ncol(gamma)), nrow=4)
  x <- c(mu_h)
  
  if(length(x) > 30) {
    q_indices <- round(length(x)/30 * seq(0, 29))
  } else {
    q_indices <- seq(0, 29)
  }
  q_indices <- q_indices+1
  q_indices <- c(q_indices, length(x))
  
  sort_i <- order(x)
  x_sorted <- x[sort_i]


  y <- 1-gamma
  y_sorted <- y[sort_i,]
  
  x_quant <- rep(0, length(q_indices)-1);
  y_quant <- matrix(0L, nrow=length(q_indices)-1, ncol=ncol(y))
  
  for(i in 1:(length(q_indices)-1)){
    start_i <- q_indices[i];
    end_i <- q_indices[i+1];
    
    x_quant[i] <- mean(x_sorted[start_i:end_i]);
    y_quant[i,] = colMeans(y_sorted[start_i:end_i,])
    print(colMeans(y_sorted[start_i:end_i,]))
  }
  
  initialGuesses <- list(c(3.5, 1), c(5.5, 1), c(1.5, .5), c(5.5, .5), c(3.5, 1.7))
  bounds <- list(c(0, Inf), c(0, 2))
  print(y_quant[,1])
  return()
  
  for (k in 1:ncol(gamma)) {
    best_eval <- 1e99
    for (ig in initialGuesses) {
      print(y_quant[,i])
      res <- optim(par=ig, fn=efun, x=x_quant, y=y_quant[,i], method="L-BFGS-B",lower = bounds[[1]], upper=bounds[[2]])
      print(res)
      return()
    }
    return("here")
      
  }
  
  
}