context("Signature Scores")

test_that(".gmt files are read in correctly", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/c2.all.v5.2.symbols.gmt"))
  expect_equal(length(sigList), 4729)
  
})

test_that(".txt files are read in correctly", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/tcga_sigs.txt"))
  expect_equal(length(sigList), 4)
  
})

test_that("Read in .txt and .gmt files together", {
  
  sigList <- readSignaturesInput(c("test_data/published_signatures/tcga_sigs.txt", 
                                    "test_data/published_signatures/c2.all.v5.2.symbols.gmt"))
  expect_equal(length(sigList), 4733)
  
})

test_that("Precomputed files read in correctly", {
  
  expr <- readExprToMatrix("test_data/expression_matrix.txt")
  prList <- readPrecomputed("test_data/precomputed_signatures.txt", colnames(expr))
  
  expect_equal(length(prList), 1)
  
  pr <- prList[[1]]
  expect_true(pr@isFactor)
  expect_true(pr@isPrecomputed)
  expect_equal(pr@sample_labels, colnames(expr))
  
})

test_that("Naive Sig Scores computed correctly", {
  
  expr <- readExprToMatrix("test_data/expression_matrix.txt")
  eData <- ExpressionData(expr)
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  
  weights <- matrix(1L, nrow=nrow(expr), ncol=ncol(expr))
  rownames(weights) <- rownames(expr)
  colnames(weights) <- colnames(expr)
  min_signature_genes <- 0
  
  sigScores <- c()
  
  for (sig in sigList) {
    tryCatch({
      sigScores <- c(sigScores, naiveEvalSignature(eData, 
                                                   sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  expect_equal(length(sigScores), length(sigList))
  ss <- sigScores[[1]]
  expect_false(ss@isPrecomputed)
  expect_false(ss@isFactor)
  expect_equal(ss@sample_labels, colnames(expr))
  expect_equal(length(ss@scores), ncol(expr))
})

test_that("Naive Sig Eval is same as Weighted Sig Eval with all weights 1", {
  expr <- readExprToMatrix("test_data/expression_matrix.txt")
  eData <- ExpressionData(expr)
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  
  weights <- matrix(1L, nrow=nrow(expr), ncol=ncol(expr))
  rownames(weights) <- rownames(expr)
  colnames(weights) <- colnames(expr)
  min_signature_genes <- 0
  
  sigScores <- c()
  
  for (sig in sigList) {
    tryCatch({
      sigScores <- c(sigScores, naiveEvalSignature(eData, 
                                                   sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  expect_equal(length(sigScores), length(sigList))
  
  sigScores2 <- c()
  for (sig in sigList) {
    tryCatch({
      sigScores2 <- c(sigScores2, weightedEvalSignature(eData, 
                                                   sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  expect_equal(length(sigScores2), length(sigList))
  
  for (i in 1:length(sigList)) {
    eq <- (sigScores[[i]]@scores == sigScores2[[i]]@scores)
    
    expect_false(FALSE %in% eq)
  }
  
})


test_that("Test Signature Scores against correct Signature Score Matrix", {
  
  t_sigMat <- as.matrix(read.table("test_data/test_sigMatrix.txt", sep="\t"))
  
  expr <- readExprToMatrix("test_data/expression_matrix.txt")
  eData <- ExpressionData(expr)
  sigList <- readSignaturesInput(c("test_data/published_signatures/h.all.v5.2.symbols.gmt"))
  housekeeping_data <- readHKGToMatrix("test_data/housekeeping_genes/Gene Name Housekeeping.txt")
  
  
  fneg_out <- createFalseNegativeMap(expr, housekeeping_data, FALSE)
  func <- fneg_out[[1]]
  params <- fneg_out[[2]]
  
  weights <- computeWeights(func, params, eData)
  
  normalizedData <- getNormalizedCopy(eData, "znorm_rows")
  eData <- updateExprData(eData, normalizedData)
  
  min_signature_genes <- 1
  sigScores <- c()
  for (sig in sigList) {
    tryCatch({
      sigScores <- c(sigScores, weightedEvalSignature(eData, 
                                                        sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  names <- c()
  sigMatrix <- matrix(0L, nrow=length(sigScores), ncol=length(sigScores[[1]]@scores))
  for (sig in 1:length(sigScores)) {
    names <- c(names, sigScores[[sig]]@name)
    scores <- t(as.matrix(sigScores[[sig]]@scores))
    sigMatrix[sig,] <- scores
  }
  
  rownames(sigMatrix) <- names
  colnames(sigMatrix) <- colnames(expr)
  
  expect_equal(colnames(sigMatrix), colnames(t_sigMat))
  expect_equal(rownames(sigMatrix), rownames(t_sigMat))
  expect_equal(sigMatrix, t_sigMat)
  
})
