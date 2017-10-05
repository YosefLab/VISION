context("Project Expression Data")

test_that("tSNE is accurate from PCA data", {

    res <- as.matrix(read.table("test_data/pca_res.txt", sep="\t"))

    tsne30 <- applytSNE30(res)
    t_tsne30 <- as.matrix(read.table("test_data/tsne30_res.txt", sep="\t"))
    colnames(t_tsne30) <- NULL

    expect_equal(tsne30, t_tsne30, tolerance=5e1)
    
})

test_that("KNN returns the same weights as computed in sigs vs projections", {

    res <- as.matrix(read.table("test_data/pca_res.txt", sep="\t"))
    proj <- Projection("pca", pData=res)

    knn <- applyKNN(res)

    knn2 <- computeKNNWeights(proj, round(sqrt(ncol(res))), SerialParam())

    expect_equal(knn, knn2)

})

