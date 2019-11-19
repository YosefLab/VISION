context("Signature Scores")

test_that(".gmt files are read in correctly", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  expect_equal(length(sigList), 48)

  # Test signature file with explicit values
  sigList <- readSignaturesInput("test_data/sigs_w_values.gmt")
  expect_equal(length(sigList), 2)

  expect_equal(sigList$sig1@sigDict, c(GENE1=0.5, GENE2=1.5))
  expect_equal(sigList$sig2@sigDict, c(GENE2=-0.5, GENE3=4))
  
  
})

test_that(".txt files are read in correctly", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/tcga_sigs.txt"))
  expect_equal(length(sigList), 4)
  
})

test_that("Read in .txt and .gmt files together", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/tcga_sigs.txt", 
                                    "test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  expect_equal(length(sigList), 52)
  
})


test_that("Naive Sig Scores computed correctly", {
  
  expr <- read.table("test_data/expression_matrix.txt",
                     sep = "\t", header = TRUE)
  expr <- data.matrix(expr)
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  sigList <- processSignatures(sigList, rownames(expr), 1)
  
  normData <- getNormalizedCopySparse(expr, "none")
  
  sigScores <- batchSigEvalNorm(sigList, normData)
  
  expect_is(sigScores, "matrix")
  expect_equal(ncol(sigScores), length(sigList))
  expect_equal(nrow(sigScores), ncol(expr))
})


test_that("Unweighted signature scores are correct", {
    expr <- c(
              1, 2.0,
              2, 2,
              5, 4
              )
    expr <- matrix(expr, nrow = 3, byrow = TRUE)
    rownames(expr) <- c("A", "B", "C")
    colnames(expr) <- c("s1", "s2")

    normData <- getNormalizedCopySparse(expr, "none")
    normData@data <- 2**normData@data - 1

    sig1 <- Signature(sigDict = c("a" = 1, "b" = 1), name = "sig1",
                      source = "", metaData = "")

    sig2 <- Signature(sigDict = c("b" = 1, "c" = -1), name = "sig2",
                      source = "", metaData = "")

    sigs <- list(sig1 = sig1, sig2 = sig2)

    scores <- batchSigEvalNorm(sigs, normData)
    expect_equal(ncol(scores), 2)

    sig1scores <- scores[, "sig1"]
    sig2scores <- scores[, "sig2"]

    expected_sig1_scores <- c(s1 = 1.5, s2 = 2.0)
    expect_equal(sig1scores, expected_sig1_scores)

    expected_sig2_scores <- c(s1 = -3.0, s2 = -2.0) / 2
    expect_equal(sig2scores, expected_sig2_scores)

})
