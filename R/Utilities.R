#' Utility methods


#' Helper utility to group list items into batches
#'
#' This is used inside the batchSigEval function
#'
#' @param items list of items to group into batches
#' @param per_batch approximate target for size of batches
#' @param n_workers number of batches is made divisible by n_workers
#' @return list of list of items
batchify <- function(items, per_batch, n_workers = 1) {
    if (length(items) == 0){
        return(list(list()))
    }
    n_iterations <- round(length(items) / (per_batch * n_workers))

    if (n_iterations == 0){
        n_batches <- ceiling(length(items) / per_batch)
    } else {
        n_batches <- n_iterations * n_workers
    }
    per_batch <- ceiling(length(items) / n_batches)

    out <- lapply(seq(n_batches), function(i) {
        start_i <- (i - 1) * per_batch + 1
        end_i <- i * per_batch
        if (end_i > length(items)){
            end_i  <- length(items)
        }
        return(items[start_i:end_i])
    })

    return(out)
}


#' Check's the version of the FastProject object and displays error if necessary
#'
#' @param object FastProject object
#' @return NULL
versionCheck <- function(object) {

    templateStr <- paste0(
        "This FastProject object was created with an older version of the library.",
        "  To view, either install an older version (commit #COMMIT) from GitHub (www.github.com/YosefLab/FastProjectR) or",
        " recreate the object and re-run analyze"
    )

    if (!.hasSlot(object, "version")) {
        msg <- gsub("#COMMIT", "0cd5268", templateStr)
        stop(msg, call. = FALSE)
    }

    if(object@version < 1.0) {
        msg <- gsub("#COMMIT", "0cd5268", templateStr)
        stop(msg, call. = FALSE)
    }

    # Add new commit hashes here as version increases and breaks backward compatibility

    return()
}


#' log2-scale transform a dense OR sparse matrix
#'
#' This avoids the creation of a dense intermediate matrix when
#' operating on sparse matrices
#'
#' Either performs result <- log2(spmat+1) or if scale = TRUE
#' returns result <- log2(spmat/colSums(spmat)*scaleFactor + 1)
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix summary
#' @param spmat sparse Matrix
#' @param scale boolean - whether or not to scale the columns to sum to `scale_factor`
#' @param scaleFactor if scale = TRUE, columns are scaled to sum to this number
#' @return logmat sparse Matrix
matLog2 <- function(spmat, scale = FALSE, scaleFactor = 1e6) {


    if (scale == TRUE) {
        spmat <- t( t(spmat) / colSums(spmat)) * scaleFactor
    }

    if (is(spmat, "sparseMatrix")) {
        matsum <- summary(spmat)

        logx <- log2(matsum$x + 1)

        logmat <- sparseMatrix(i = matsum$i, j = matsum$j,
                               x = logx, dims = dim(spmat),
                               dimnames = dimnames(spmat))
    } else {
        logmat <- log2(spmat + 1)
    }


    return(logmat)

}

#' Vectorized wilcox rank-sums test
#'
#' Given indices in cluster_ii, compute the ranksums test
#' on every column for values in cluster_ii vs the rest.
#'
#' @importFrom matrixStats colCounts
#' @importFrom matrixStats colMaxs
#' @param ranks matrix of ranks, each column representing a separate variable
#' @param cluster_ii numeric vector - indices denoting the group to be compared
#' @param check_na - whether or not to check for NA values (slows computation)
#' @param check_ties - whether or not to check for ties in ranks (slows computation)
#' @return pval - numeric vector, pvalue for each column
#' @return stat - numeric vector, test statistic (scaled AUC) for each column
matrix_wilcox <- function(ranks, cluster_ii,
                          check_na = FALSE, check_ties = FALSE){

    # handle edge case
    if (ncol(ranks) == 0){
        p <- numeric()
        stat <- numeric()
        return(list(pval = p, stat = stat))
    }

    subset <- ranks[cluster_ii, , drop = FALSE]

    not_cluster_ii = setdiff(seq(nrow(ranks)), cluster_ii)
    subset_2 = ranks[not_cluster_ii, , drop=FALSE]

    if (check_na){
        n1 <- nrow(subset) - colCounts(subset, value = NA)
	n2 <- nrow(subset_2) - colCounts(subset_2, value=NA)
    } else {
        n1 <- length(cluster_ii)
        n2 <- nrow(ranks) - n1
    }


    r1 <- colSums(subset, na.rm = T)

    u1 <- r1 - n1 * (n1 + 1) / 2

    u <- u1
    AUC <- u / (n1 * n2)
    AUC <- ifelse(is.infinite(AUC), .5, AUC) # Edge case n2 == 0
    AUC <- ifelse(is.na(AUC), .5, AUC) # Edge case, n1 == 0
    stat <- (AUC - .5) * 2  # translates (0.5, 1) to (0, 1)

    # u to z-score
    m_u <- (n1 * n2) / 2
    if (check_ties){
        sd_u <- vapply(seq_len(ncol(ranks)), function(i){
                    Y <- ranks[, i]
                    n1 <- if (length(n1) == 1) n1 else n1[i]
                    n2 <- if (length(n2) == 1) n2 else n2[i]
                    has_ties <- length(Y) != length(unique(Y))
                    if (!has_ties){
                        sd_ui <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
                        return(sd_ui)
                    }
                    n <- n1 + n2
                    n_ties <- table(Y)
                    sd_ui <- sqrt(
                                  n1 * n2 / 12 *
                                  (n + 1 -
                                   sum(n_ties ** 3 - n_ties) / n / (n - 1)
                                  )
                             )
                    return(sd_ui)
        }, FUN.VALUE = 0.0)


    } else {
        sd_u <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    }
    z <- (u - m_u) / sd_u
    z <- -1 * abs(z) # ensure negative
    z <- z + .5 / sd_u # continuity correction
    z[sd_u == 0] <- 0  # handle case where n1 or n2 = 0
    p <- pnorm(z) * 2

    p[p > 1] <- 1  # Can happen due to continuity correction

    return(list(pval = p, stat = stat))
}

#' Perform 1vAll factor analysis given a factor matrix and group definition
#'
#' Given indices in cluster_ii, compute the chisq test
#' on every column for values in cluster_ii vs the rest.
#'
#' @param factorDF dataframe of factors
#' @param cluster_ii numeric vector - indices denoting the group to be compared
#' @return pval - numeric vector, pvalue for each column
#' @return stat - numeric vector, test statistic (Cramer's V) for each column
matrix_chisq <- function(factorDF, cluster_ii) {

    out <- lapply(colnames(factorDF), function(var){

                    if (!is.factor(factorDF[[var]])){
                        stop("Error: matrix_chisq must be called on a factor dataframe")
                    }

                    values <- factorDF[, var, drop = F]

                    values[, 2] <- 0
                    values[cluster_ii, 2] <- 1
                    values[, 2] <- as.factor(values[, 2])
                    M <- table(values)
                    if (nrow(M) > 1 && ncol(M) > 1){
                        suppressWarnings(
                            out <- chisq.test(M)
                        )
                        pval <- out$p.value
                        n <- sum(M)
                        V <- sqrt(out$statistic / n /
                            min(nrow(M) - 1, ncol(M) - 1)
                        )
                        stat <- V # Cramer's V
                    } else {
                        pval <- 1.0
                        stat <- 0
                    }
                return(list(pval = pval, stat = stat))
    })
    names(out) <- colnames(factorDF)

    pvals <- vapply(out, function(x) x$pval, FUN.VALUE = 0.0)
    stat <- vapply(out, function(x) x$stat, FUN.VALUE = 0.0)

    return(list(pval = pvals, stat = stat))
}
