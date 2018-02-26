context("Project Expression Data")

test_that("KNN returns the same weights as computed in sigs vs projections", {

    expr <- readExprAsMatrix("test_data/expression_matrix.txt")
    res <- applyPCA(expr, N=10)[[1]]
    proj <- Projection("pca", pData=res)

    knn <- applyKNN(res)

    knn2 <- computeKNNWeights(proj, round(sqrt(ncol(res))))

    total_diff <- sum(abs(knn - knn2))
    expect_equal(total_diff, 0)

})

