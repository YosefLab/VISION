context("Creating Weights & Background Distribution")

test_that("Testing weight matrix created correctly", {
  
  expr <- readExprToMatrix("test_data/expression_matrix.txt")
  eData <- ExpressionData(expr)
  housekeeping_data <- readHKGToMatrix("test_data/housekeeping_genes/Gene Name Housekeeping.txt")
  
  
  fneg_out <- createFalseNegativeMap(expr, housekeeping_data, FALSE)
  func <- fneg_out[[1]]
  params <- fneg_out[[2]]
  
  expect_equal(ncol(params), ncol(expr))
  expect_equal(nrow(params), 4)
  
  weights <- computeWeights(func, params, eData)
  
  expect_equal(ncol(weights), ncol(expr))
  expect_equal(nrow(weights), nrow(expr))
  
  # Test that all values are between 0 and 1
  t_weights <- apply(weights, 1, function(x) (x >= 0 && x <= 1))
  expect_false(FALSE %in% t_weights)
})

test_that("Compare Output to Previous Weights Matrix", {
  expr <- readExprToMatrix("test_data/expression_matrix.txt")
  eData <- ExpressionData(expr)
  housekeeping_data <- readHKGToMatrix("test_data/housekeeping_genes/Gene Name Housekeeping.txt")
  
  
  fneg_out <- createFalseNegativeMap(expr, housekeeping_data, FALSE)
  func <- fneg_out[[1]]
  params <- fneg_out[[2]]
  
  weights <- computeWeights(func, params, eData)
  t_weights <- as.matrix(read.table("test_data/test_weights.txt", sep="\t"))
  
  expect_equal(weights, t_weights)
  
  
})