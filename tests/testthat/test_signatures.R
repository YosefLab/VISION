context("Signature Scores")

test_that(".gmt files are read in correctly", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  expect_equal(length(sigList), 48)
  
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
  
  expr <- readExprAsMatrix("test_data/expression_matrix.txt")
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  
  weights <- matrix(1L, nrow=nrow(expr), ncol=ncol(expr))
  rownames(weights) <- rownames(expr)
  colnames(weights) <- colnames(expr)
  
  sigScores <- batchSigEval(sigList, "naive", expr, weights, 1)
  
  expect_equal(length(sigScores), length(sigList)-1)
  ss <- sigScores[[1]]
  expect_false(ss@isMeta)
  expect_false(ss@isFactor)
  expect_equal(names(ss@scores), colnames(expr))
  expect_equal(length(ss@scores), ncol(expr))
})

test_that("Naive Sig Eval is same as Weighted Sig Eval with all weights 1", {
  expr <- readExprAsMatrix("test_data/expression_matrix.txt")
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  
  weights <- matrix(1L, nrow=nrow(expr), ncol=ncol(expr))
  rownames(weights) <- rownames(expr)
  colnames(weights) <- colnames(expr)

  sigScores <- batchSigEval(sigList, "naive", expr, weights, 1)
  
  expect_equal(length(sigScores), length(sigList)-1)
  
  sigScores2 <- batchSigEval(sigList, "weighted_avg", expr, weights, 1)
  
  expect_equal(length(sigScores2), length(sigList)-1)
  
  for (i in 1:length(sigScores2)) {
    eq <- (sigScores[[i]]@scores == sigScores2[[i]]@scores)
    
    expect_false(FALSE %in% eq)
  }
  
})

test_that("Weighted signature scores are correct", {
    expr <- c(
              1, 2.0,
              2, 2,
              5, 4
              )
    expr <- matrix(expr, nrow = 3, byrow = TRUE)
    rownames(expr) <- c("A", "B", "C")
    colnames(expr) <- c("s1", "s2")

    weights <- c(
              .5, .3,
              .1, .5,
              .3, .6
              )
    weights <- matrix(weights, nrow = 3, byrow = TRUE)
    rownames(weights) <- rownames(expr)
    colnames(weights) <- colnames(expr)

    sig1 <- Signature(sigDict = c("a" = 1, "b" = 1), name = "sig1",
                      source = "", metaData = "", isMeta = FALSE,
                      isFactor = FALSE)

    sig2 <- Signature(sigDict = c("b" = 1, "c" = -1), name = "sig2",
                      source = "", metaData = "", isMeta = FALSE,
                      isFactor = FALSE)

    sigs <- list(sig1 = sig1, sig2 = sig2)

    scores <- batchSigEval(sigs, "weighted_avg", expr, weights, 1)
    expect_equal(length(scores), 2)

    sig1scores <- scores[["sig1"]]


    sig2scores <- scores[["sig2"]]

    expected_sig1_scores <- c(s1 = (1 * .5 + 2 * .1) / (.5 + .1),
                              s2 = (2 * .3 + 2 * .5) / (.3 + .5))
    expect_equal(sig1scores@scores, expected_sig1_scores)
    expect_equal(sig1scores@name, sig1@name)

    expected_sig2_scores <- c(s1 = (2 * .1 - 5 * .3) / (.1 + .3),
                              s2 = (2 * .5 - 4 * .6) / (.5 + .6))
    expect_equal(sig2scores@scores, expected_sig2_scores)
    expect_equal(sig2scores@name, sig2@name)

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

    weights <- matrix(1.0, nrow = nrow(expr), ncol = ncol(expr))
    rownames(weights) <- rownames(expr)
    colnames(weights) <- colnames(expr)

    sig1 <- Signature(sigDict = c("a" = 1, "b" = 1), name = "sig1",
                      source = "", metaData = "", isMeta = FALSE,
                      isFactor = FALSE)

    sig2 <- Signature(sigDict = c("b" = 1, "c" = -1), name = "sig2",
                      source = "", metaData = "", isMeta = FALSE,
                      isFactor = FALSE)

    sigs <- list(sig1 = sig1, sig2 = sig2)

    scores <- batchSigEval(sigs, "naive", expr, weights, 1)
    expect_equal(length(scores), 2)

    sig1scores <- scores[["sig1"]]


    sig2scores <- scores[["sig2"]]

    expected_sig1_scores <- c(s1 = 1.5, s2 = 2.0)
    expect_equal(sig1scores@scores, expected_sig1_scores)
    expect_equal(sig1scores@name, sig1@name)

    expected_sig2_scores <- c(s1 = -3.0, s2 = -2.0) / 2
    expect_equal(sig2scores@scores, expected_sig2_scores)
    expect_equal(sig2scores@name, sig2@name)

})
