context("Matrix Wilcox Test")

test_that("Calculation is correct, no NA", {

    test_mat <- matrix(runif(5*10), nrow=5)
    test_mat_rank <- colRanks(test_mat,
                              ties.method = "average",
                              preserveShape = TRUE)
    cluster_ii <- c(1, 2, 5)
    not_cluster_ii <- which(!(seq(nrow(test_mat)) %in% cluster_ii))

    # Our implementation
    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = FALSE, check_ties = TRUE)
    p_test <- out$pval
    stat_test <- out$stat

    # Compare vs R's
    out <- lapply(seq_len(ncol(test_mat)), function(i){
                      a <- test_mat[cluster_ii, i]
                      b <- test_mat[not_cluster_ii, i]
                      n1 <- sum(!is.na(a))
                      n2 <- sum(!is.na(b))
                      test <- wilcox.test(a, b, alternative = "two.sided", exact = FALSE)
                      AUC <- test$statistic / (n1 * n2)
                      return(list(pval = test$p.value, stat = AUC))
    })

    p_correct <- vapply(out, function(x) x$pval, FUN.VALUE = 0.0)
    stat_correct <- vapply(out, function(x) x$stat, FUN.VALUE = 0.0)

    expect_equal(p_test, p_correct, tolerance = 1e-8)
    expect_equal(stat_test, stat_correct, tolerance = 1e-8)
    expect_equal(length(p_test), ncol(test_mat))
    expect_equal(length(stat_test), ncol(test_mat))
})

test_that("Calculation is correct with ties, no NA", {

    test_mat <- matrix(runif(5*10), nrow=5)
    test_mat[c(2, 3, 4), 5] <- 1 # Add some ties
    test_mat_rank <- colRanks(test_mat,
                              ties.method = "average",
                              preserveShape = TRUE)
    cluster_ii <- c(1, 2, 5)
    not_cluster_ii <- which(!(seq(nrow(test_mat)) %in% cluster_ii))

    # Our implementation
    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = FALSE, check_ties = TRUE)
    p_test <- out$pval
    stat_test <- out$stat

    # Compare vs R's
    out <- lapply(seq_len(ncol(test_mat)), function(i){
                      a <- test_mat[cluster_ii, i]
                      b <- test_mat[not_cluster_ii, i]
                      n1 <- sum(!is.na(a))
                      n2 <- sum(!is.na(b))
                      test <- wilcox.test(a, b, alternative = "two.sided", exact = FALSE)
                      AUC <- test$statistic / (n1 * n2)
                      return(list(pval = test$p.value, stat = AUC))
    })

    p_correct <- vapply(out, function(x) x$pval, FUN.VALUE = 0.0)
    stat_correct <- vapply(out, function(x) x$stat, FUN.VALUE = 0.0)

    expect_equal(p_test, p_correct, tolerance = 1e-8)
    expect_equal(stat_test, stat_correct, tolerance = 1e-8)
})

test_that("Calculation is correct with NAs", {

    test_mat <- matrix(runif(5*10), nrow=5)
    # Add some NAs
    test_mat[1, 1] <- NA
    test_mat[1, 2] <- NA
    test_mat[3, 2] <- NA
    test_mat[3, 3] <- NA
    test_mat_rank <- colRanks(test_mat,
                              ties.method = "average",
                              preserveShape = TRUE)
    cluster_ii <- c(1, 2, 5)
    not_cluster_ii <- which(!(seq(nrow(test_mat)) %in% cluster_ii))

    # Our implementation
    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = TRUE, check_ties = TRUE)
    p_test <- out$pval
    stat_test <- out$stat

    # Compare vs R's
    out <- lapply(seq_len(ncol(test_mat)), function(i){
                      a <- test_mat[cluster_ii, i]
                      b <- test_mat[not_cluster_ii, i]
                      n1 <- sum(!is.na(a))
                      n2 <- sum(!is.na(b))
                      test <- wilcox.test(a, b, alternative = "two.sided", exact = FALSE)
                      AUC <- test$statistic / (n1 * n2)
                      return(list(pval = test$p.value, stat = AUC))
    })

    p_correct <- vapply(out, function(x) x$pval, FUN.VALUE = 0.0)
    stat_correct <- vapply(out, function(x) x$stat, FUN.VALUE = 0.0)

    expect_equal(p_test, p_correct, tolerance = 1e-8)
    expect_equal(stat_test, stat_correct, tolerance = 1e-8)
    expect_equal(length(p_test), ncol(test_mat))
    expect_equal(length(stat_test), ncol(test_mat))
})

test_that("Edge cases produce expected outputs", {

    #   single col
    single_col <- matrix(seq_len(5), ncol = 1)

    cluster_ii <- c(1, 2, 3)

    out <- matrix_wilcox(single_col, cluster_ii,
                         check_na = FALSE, check_ties = FALSE)

    expect_equal(length(out$pval), 1)
    expect_equal(length(out$stat), 1)

    out <- matrix_wilcox(single_col, cluster_ii,
                         check_na = TRUE, check_ties = TRUE)

    expect_equal(length(out$pval), 1)
    expect_equal(length(out$stat), 1)

    # zero columns
    zero_col <- matrix(nrow = 5, ncol = 0)

    out <- matrix_wilcox(zero_col, cluster_ii,
                         check_na = FALSE, check_ties = FALSE)

    expect_equal(length(out$pval), 0)
    expect_equal(length(out$stat), 0)

    out <- matrix_wilcox(zero_col, cluster_ii,
                         check_na = TRUE, check_ties = TRUE)

    expect_equal(length(out$pval), 0)
    expect_equal(length(out$stat), 0)

    # empty cluster_ii

    test_mat <- matrix(runif(5*10), nrow=5)
    test_mat_rank <- colRanks(test_mat,
                              ties.method = "average",
                              preserveShape = TRUE)
    cluster_ii <- c()

    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = FALSE, check_ties = FALSE)

    expect_equal(out$pval, rep(1, 10))
    expect_equal(out$stat, rep(0.5, 10))

    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = TRUE, check_ties = TRUE)

    expect_equal(out$pval, rep(1, 10))
    expect_equal(out$stat, rep(0.5, 10))

    # full cluster_ii
    cluster_ii <- seq(5)

    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = FALSE, check_ties = FALSE)

    expect_equal(out$pval, rep(1, 10))
    expect_equal(out$stat, rep(0.5, 10))

    out <- matrix_wilcox(test_mat_rank, cluster_ii,
                         check_na = TRUE, check_ties = TRUE)

    expect_equal(out$pval, rep(1, 10))
    expect_equal(out$stat, rep(0.5, 10))

})

test_that("CPP Wilcox is Correct", {

    test_mat <- matrix(runif(5 * 10), nrow = 5)
    test_mat[3, c(3, 4, 7)] <- 1 # Add some ties
    test_mat[4, c(1, 2, 4)] <- 1 # Add some ties
    test_mat[5, c(2, 3, 4)] <- 1 # Add some ties
    cluster_ii <- c(1, 2, 5)
    not_cluster_ii <- which(!(seq(ncol(test_mat)) %in% cluster_ii))

    res <- matrix_wilcox_cpp(test_mat, cluster_ii, not_cluster_ii)


    # Compare with R implementation
    ref <- do.call(rbind,
        lapply(seq_len(nrow(test_mat)), function(i){
            x <- test_mat[i, cluster_ii]
            y <- test_mat[i, not_cluster_ii]
            ref_i <- wilcox.test(x, y, exact = FALSE, alternative = "two.sided")
            return(c(U = unname(ref_i$statistic), pval = ref_i$p.value))
        })
    )

    ref <- as.data.frame(ref)

    expect_equal(ref$U, res$U)
    expect_equal(ref$pval, res$pval, tolerance = 1e-8)
})
