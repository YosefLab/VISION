## Functions to transform the data and calculate weights
require(loe)

louvainCluster <- function(kn, data) {
	
	nn <- kn[[1]]
	d <- kn[[2]]
	nnl <- lapply(1:nrow(nn), function(i) nn[i,])

	# Create an undirected knn graph
	g <- graph_from_adj_list(nnl, mode="out")
	E(g)$weights <- as.vector(t(d))
	g <- as.undirected(g, mode="each")
	
	# Now apply the louvain algorithm to cluster the graph
	cl <- cluster_louvain(g)

	# Gather cluster vector to list of clusters
	clusters <- list()
	mem <- as.vector(membership(cl))
	for (i in 1:length(mem)) {
		n <- as.character(mem[[i]])
		if (n %in% names(clusters)) {
			clusters[[n]] <- c(clusters[[n]], i)
		} else {
			clusters[[n]] <- c(i)
		}
	}

	clusters <- lapply(clusters, function(i) i <- rownames(data)[i])


	return(clusters)

}

readjust_clusters <- function(clusters, data) {
	#' Repartitions existing clusters to achieve desired granularity
	#' Paramters:
	#'	clusters: (List) list of clusters, each entry begin a vector of cells in a cluster
	#'	data: (matrix) NUM_SAMPLES x NUM_FEATURES data matrix that was used to generate clusters

	NUM_PARTITIONS = sqrt(nrow(data))
	EPSILON = .15

	currPart = length(clusters)
	clusterList <- list()
	
	while (currPart < ((1 - EPSILON)*NUM_PARTITIONS)) {
		clusterList <- list()
		cluster_offset = 0
		for (i in 1:length(clusters)) {

			# Apply kmeans clustering to existing cluster
			currCl = clusters[[i]]
			subData <- data[currCl,]
			nCl <- kmeans(subData, centers=round(sqrt(length(currCl))), iter.max=100)
			
			# Gather cluster vector to list of clusters
			for (i in 1:length(nCl$cluster)) {
				n <- as.character(nCl$cluster[[i]] + cluster_offset)
				sample_n <- names(nCl$cluster)[[i]]
				if (n %in% names(clusterList)) {
					clusterList[[n]] <- c(clusterList[[n]], sample_n)
				} else {
					clusterList[[n]] <- c(sample_n)
				}

			}

			# Now add to cluster offset for next re-clustering
			cluster_offset <- cluster_offset + max(nCl$cluster)
		}

		currPart <- length(clusterList)
		clusters <- clusterList

	}

	return(clusters)
}


createFalseNegativeMap <- function(data, housekeeping_genes, debug=0) {
  #' Uses gene names in `housekeeping_genes` to create a mapping of false negatives.
  #' Creates a functional fit for each sample based on that samples HK genes
  #'
  #' Returns: fit_func: (function) Used to fit expression values to FN rate
  #'          params: (data.frame - Num_Params x Num_Samples) 
  #'                  Sample-specific parameters to use with fit_func
 
  message("Creating False Negative Map...")
  
  # get subset of genes to be used, ie those included in the housekeeping genes set
  keep_ii <- which(rownames(data) %in% housekeeping_genes[[1]])
  
  # Filter out genes with no variance
  data_hk <- data[keep_ii, ]
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

 
  params <- matrix(0L, ncol=ncol(gamma), nrow=4)
  
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
        param <- res$par
        params[1,k] = param[1]
        params[2,k] = param[2]
        params[3,k] = 0
        params[4,k] = 1
        
      }
    }
  }
  
  if (debug > 0) {
    i = debug;
    domain=seq(0,10,0.001)
    
    plot(domain, func(domain, params[1,i], params[2,i], params[3,i], params[4,i]), 
        type='l', xlim=c(0, 10), ylim=c(0, 1), 
        ylab = paste0("P(gene not expressed in ", colnames(data)[i], ")"),
        xlab = "Gene average in samples expressing gene")
    points(x,y[,i], 'p', col='blue');
    points(x_quant, y_quant[,i], 'p', col='red')
    print(params[,i])
  }
  
  return(list(func, params))
  
}

computeWeights <- function(fit_func, params, exprData) {
  #' Calculates weights for the data from the FNR curves
  #' Weights represent p(not expressed | not detected) for zero values
  #' and are equal to 1.0 for detected values.
  #' Parameters: fit_func (function) (mu_h, params) 
  #'             Function parametrized by params that maps each mu_h to a false negative estimate
  #'             
  #'             params (4 x Num_Samples Matrix) 
  #'             Matrix containing parameters for the false negative  fit function (fit_func)
  #'             
  #'             expr (ExpressionData) 
  #'             Data from which prob derives
  #'             
  #' Returns: weights (Num_Genes x Num_Samples Matrix) Estimated weight for each data point in input matrix
  #'                  Ranges from 0 to 1. 
  
  expr <- getExprData(exprData);

  fnProb <- matrix(0L, nrow = nrow(expr), ncol = ncol(expr))
  countNonZero <- apply(expr, 1, function(c) sum(c!=0))
  countNonZero[which(countNonZero == 0)] <- 1;
  mu_h <- apply(expr, 1, function(r) sum(r)) / countNonZero
  
  for (i in 1:ncol(fnProb)) {
    fnProb[,i] = fit_func(mu_h, params[,i][1], params[,i][2], params[,i][3], params[,i][4])
  }
  
  pdE <- 1 - fnProb
  pnd <- apply(expr, 1, function(r) sum(r==0)) / ncol(expr)
  pe <- (1 - pnd) / apply(pdE, 1, function(r) mean(r))
  
  pe[which(is.na(pe))] <- 1.0
  pnd[which(pnd == 0)] <- 1.0 / ncol(expr)
  
  pne_nd <- 1 - (1-pdE)* (pe / pnd)
  pne_nd[which(pne_nd < 0)] <- 0.0
  pne_nd[which(pne_nd > 1)] <- 1.0
  
  weights <- pne_nd
  weights[which(expr > 0)] <- 1.0

  rownames(weights) <- rownames(expr)
  colnames(weights) <- colnames(expr)
  
  return(weights)
  
}
