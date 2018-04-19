context("Geary C Calculation")

test_that("Calculation is correct on dense", {
    X <- c(1, 1, 1, 1, 2, 2, 2, 2)
    W <- c(
           1, 1, 1, 1, .1, .1, .1, .1,
           1, 1, 1, 1, .1, .1, .1, .1,
           1, 1, 1, 1, .1, .1, .1, .1,
           1, 1, 1, 1, .1, .1, .1, .1,
           .1, .1, .1, .1, 1, 1, 1, 1,
           .1, .1, .1, .1, 1, 1, 1, 1,
           .1, .1, .1, .1, 1, 1, 1, 1,
           .1, .1, .1, .1, 1, 1, 1, 1
           )
    W <- matrix(W, nrow = 8, ncol = 8)

    result1 <- geary(X, W)

    expect_equal(result1, 0.1590909, tolerance = 1e-6)

    X2 <- c(1, 2, 1, 2, 1, 2, 1, 2)
    result2 <- geary(X2, W)
    expect_equal(result2, 0.875, tolerance = 1e-6)

    Xmat <- matrix(c(X, X2), nrow = 2, byrow = TRUE)
    result_mat <- geary_all(Xmat, W)

    expect_equal(result1, result_mat[1])
    expect_equal(result2, result_mat[2])

})

test_that("Calculation is correct on sparse", {
    X <- c(1, 1, 1, 1, 2, 2, 2, 2)
    X2 <- c(1, 2, 1, 2, 1, 2, 1, 2)
    Xmat <- matrix(c(X, X2), nrow = 2, byrow = TRUE)

    W <- c(
           1, 1, 1, 1, 0, 0, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1
           )
    W <- matrix(W, nrow = 8, ncol = 8)

    result_dense <- geary_all(Xmat, W)


    ind <- c(
             1, 2, 3, 4,
             1, 2, 3, 4,
             1, 2, 3, 4,
             1, 2, 3, 4,
             5, 6, 7, 8,
             5, 6, 7, 8,
             5, 6, 7, 8,
             5, 6, 7, 8
            )
    ind <- matrix(ind, nrow = 8, byrow = TRUE)

    indW <- c(
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1
             )
    indW <- matrix(indW, nrow = 8, byrow = TRUE)

    result1 <- geary_sparse(X, ind, indW)
    result2 <- geary_sparse(X2, ind, indW)

    expect_equal(result1, result_dense[1])
    expect_equal(result2, result_dense[2])

    result_sparse <- geary_sparse_all(Xmat, ind, indW)

    expect_equal(result1, result_sparse[1])
    expect_equal(result2, result_sparse[2])

})
