#' Uses gene names in the housekeeping genes file to create a mapping of false
#' negatives. Creates a functional fit for each sample based on that sample's
#' HK genes
#' @importFrom stats optim
#' @param data Data matix
#' @param housekeeping_genes Housekeeping gene table
#' @return Fit function used to fit expression values to FN rate
#' @return Sample specific parameters to use with the fit function
createFalseNegativeMap <- function(data, housekeeping_genes) {

    #subset of genes to be used,ie those included in the housekeeping genes set
    data_hk <- data[rownames(data) %in% housekeeping_genes,]
    # Filter out genes with no variance
    data_hk <- data_hk[filterGenesNovar(data_hk), ]

    # calculate the distributions for hk gene
    # Gamma is 1 for any non-zero data point
    # Mu_h is the row (per gene) average of non zero points
    gamma <- as.matrix(data_hk) > 0
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
    q_indices <- c(q_indices, length(x))
    } else {
    q_indices <- seq(0, length(mu_h))
    }


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

    x_quant[i] <- mean(x_sorted[start_i:min(end_i, length(x_sorted))], );
    y_quant[i,] = colMeans(as.matrix(y_sorted[start_i:min(end_i, length(y_sorted))], ))

    }

    bounds <- list(c(0, Inf), c(0, 2))
    initialGuesses <- list(c(3.5, 1),
                            c(5.5, 1),
                            c(1.5, .5),
                            c(5.5, .5),
                            c(3.5, 1.7))

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

    return(list(func, params))

}

#' Calculates weights for the data from the FNR curves
#' Weights represent p(not expressed | not detectd) for zero values and are
#' equal to 1.0 for detected values
#'
#' @param fit_func Function parameterized by params that maps each mu_h to a
#' false negative estimate
#' @param params (4 x NUM_SAMPLES) Matrix containing parameters for the false
#' negative fit function
#' @param exprData Data from which probability derives
#' @return Weight matrix (NUM_GENES x NUM_SAMPLES) which includes the estimated
#' weight for each data point in input matrix. Ranges form 0 to 1.
computeWeights <- function(fit_func, params, exprData) {
    expr <- exprData

    fnProb <- matrix(0L, nrow = nrow(expr), ncol = ncol(expr))
    countNonZero <- apply(expr, 1, function(c) sum(c!=0))
    countNonZero[countNonZero == 0] <- 1
    mu_h <- apply(expr, 1, function(r) sum(r)) / countNonZero

    for (i in 1:ncol(fnProb)) {
    fnProb[,i] = fit_func(mu_h, params[,i][1],
                            params[,i][2],
                            params[,i][3],
                            params[,i][4])
    }

    pdE <- 1 - fnProb
    pnd <- apply(expr, 1, function(r) sum(r==0)) / ncol(expr)
    pe <- (1 - pnd) / apply(pdE, 1, function(r) mean(r))

    pe[is.na(pe)] <- 1.0
    pnd[pnd == 0] <- 1.0 / ncol(expr)

    pne_nd <- 1 - (1-pdE)* (pe / pnd)
    pne_nd[pne_nd < 0] <- 0.0
    pne_nd[pne_nd > 1] <- 1.0

    weights <- pne_nd
    weights[expr > 0] <- 1.0

    rownames(weights) <- rownames(expr)
    colnames(weights) <- colnames(expr)

    return(weights)

}
