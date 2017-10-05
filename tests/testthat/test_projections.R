context("Project Expression Data")

test_that("KNN returns the same weights as computed in sigs vs projections", {

    expr <- readExprToMatrix("test_data/expression_matrix.txt")
    res <- applyPCA(expr, N=10)[[1]]
    proj <- Projection("pca", pData=res)

    knn <- applyKNN(res)

    knn2 <- computeKNNWeights(proj, round(sqrt(ncol(res))), SerialParam())

    expect_equal(knn, knn2)

})

