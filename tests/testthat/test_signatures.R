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
  min_signature_genes <- 0
  
  sigScores <- c()
  
  for (sig in sigList) {
    tryCatch({
      sigScores <- c(sigScores, naiveEvalSignature(expr,
                                                   sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  expect_equal(length(sigScores), length(sigList)-1)
  ss <- sigScores[[1]]
  expect_false(ss@isPrecomputed)
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
  min_signature_genes <- 0
  
  sigScores <- c()
  
  for (sig in sigList) {
    tryCatch({
      sigScores <- c(sigScores, naiveEvalSignature(expr,
                                                   sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  expect_equal(length(sigScores), length(sigList)-1)
  
  sigScores2 <- c()
  for (sig in sigList) {
    tryCatch({
      sigScores2 <- c(sigScores2, weightedEvalSignature(expr,
                                                   sig, weights, min_signature_genes))
    }, error=function(e){})
  }
  
  expect_equal(length(sigScores2), length(sigList)-1)
  
  for (i in 1:length(sigScores2)) {
    eq <- (sigScores[[i]]@scores == sigScores2[[i]]@scores)
    
    expect_false(FALSE %in% eq)
  }
  
})
