## Functions to transform the data and calculate weights
require(loe)

#' Applies the Louvain algorithm to generate micro-clustered data
#' 
#' @param kn List of nearest neighbor indices and euclidean distances to these nearest neighbors
#' @param data Data matrix
#' @return List of clusters, each entry being a vector of indices representing samples in the cluster.
louvainCluster <- function(kn, data) {
	
	nn <- kn[[1]]
	d <- kn[[2]]

    d <- exp(-1 * (d*d) / .05)

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

#' Repartitions existing clusters to achieve desired granularity. By default, minimum number of clusters
#' to be generated is the squareroot of the number of cells.
#' 
#' @param clusters List of clusters, each entry being a vector of cells in a cluster.
#' @param data NUM_SAMPLES x NUM_FEATURES data matrix that was used to generate clusters
#' @return Repartitioned clusters, such that a desireable number of microclusters is acheived. 
readjust_clusters <- function(clusters, data, cellsPerPartition=100) {
	#' Repartitions existing clusters to achieve desired granularity
	#' Args:
	#'	clusters: (List) list of clusters, each entry begin a vector of cells in a cluster
	#'	data: (matrix) NUM_SAMPLES x NUM_FEATURES data matrix that was used to generate clusters
	#
	#	 Returns:
	#		Repartitioned clusters, such that a desireable number of microclusters is achieved. 

	NUM_PARTITIONS = round(nrow(data) / cellsPerPartition)
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
			if (length(currCl) > cellsPerPartition) {
				nCl <- kmeans(subData, centers=round(nrow(subData) / cellsPerPartition), iter.max=100)
			} else {
				nCl <- kmeans(subData, centers=1, iter.max=100)
			}
			newClust <- nCl$cluster
			# Gather cluster vector to list of clusters
			for (i in 1:length(newClust)) {
				n <- as.character(newClust[[i]] + cluster_offset)
				sample_n <- names(newClust)[[i]]
				if (n %in% names(clusterList)) {
					clusterList[[n]] <- c(clusterList[[n]], sample_n)
				} else {
					clusterList[[n]] <- c(sample_n)
				}

			}

			# Now add to cluster offset for next re-clustering
			cluster_offset <- cluster_offset + max(newClust)
		}

		currPart <- length(clusterList)
		clusters <- clusterList

	}

	return(clusters)
}

findRepSubset <- function(ls) {
	
	ls <- unique(sort(unlist(ls)))

	intervals <- lapply(as.list(as.numeric(ls)), function(x) return(list(x - x/5, x + x/5)))
	n_intervals <- merge_intervals(intervals)

	reprsub <- lapply(n_intervals, function(i) round((i[[1]] + i[[2]]) / 2))
	
	return(reprsub)

}

merge_intervals <- function(intervals) {
	
	for (i in 1:length(intervals)) {
		interval <- intervals[[i]]
		if (i < length(intervals) - 1) {
			nxt <- intervals[[i+1]]
			if (interval[[2]] >= nxt[[1]]) {
				n_int = list(nxt[[1]], interval[[2]])
				intervals[[i+1]] <- n_int
				intervals[[i]] <- NULL
				return(merge_intervals(intervals))		
			}
		}
	}

	return(intervals)
	
}

createPools <- function(cl, expr, hkg=list()) {

	pooled_cells <- matrix(unlist(lapply(cl, function(clust) {

			clust_data <- expr[,clust]
            if (is.null(dim(clust_data))) {
                return(clust_data)
            }   
            if (length(hkg) > 0) {
	            fnr <- createFalseNegativeMap(clust_data, hkg)
                clust_weights <- computeWeights(fnr[[1]], fnr[[2]], ExpressionData(clust_data))
			    p_cell <- as.matrix(apply(clust_data, 1, sum))
                cell_norm <- as.matrix(apply(clust_weights, 1, sum))

			    return(p_cell / cell_norm)
            }

            p_cell <- as.matrix(apply(clust_data, 1, mean))
            return(p_cell)
	  })), nrow=nrow(expr), ncol=length(cl))

    return(pooled_cells)

}

#' Uses gene names in the housekeeping genes file to create a mapping of false negatives.
#' Creates a functional fit for each sample based on that sample's HK genes
#' 
#' @param data Data matix
#' @param houskeeping_genes Housekeeping gene table
#' @param debug If 0, don't enable debugging; if 1, plot FNR. Default is 0.
#' @return Fit function used to fit expression values to FN rate
#' @return Sample specific parameters to use with the fit function
createFalseNegativeMap <- function(data, housekeeping_genes, debug=0) {

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
    if (args[[1]] < 0) {
      args[[1]] = 0
    } else if (args[[1]] > Inf) {
      args[[1]] = Inf
    }
    
    if (args[[2]] < 0) {
      args[[2]] = 0
    } else if (args[[2]] > 2) {
      args[[2]] = 2
    }
    out <- func(x, args[[1]], args[[2]])
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

#' Calculates weights for the data from the FNR curves
#' Weights represent p(not expressed | not detectd) for zero values and are equal to 1.0 for detected values
#' 
#' @param fit_func Function parameterized by params that maps each mu_h to a false negative estimate
#' @param params (4 x NUM_SAMPLES) Matrix containing parameters for the false negative fit function
#' @param exprData Data from which probability derives
#' @return Weight matrix (NUM_GENES x NUM_SAMPLES) which includes the estimated weight for each data point in
#' input matrix. Ranges form 0 to 1.
computeWeights <- function(fit_func, params, exprData) {
  #' Calculates weights for the data from the FNR curves
  #' Weights represent p(not expressed | not detected) for zero values
  #' and are equal to 1.0 for detected values.
	#'
  #' Args:			fit_func (function) (mu_h, params) 
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
