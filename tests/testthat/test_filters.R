context("Filter Expression Data")

test_that("No filter returns same matrix", {
  data <- readExprToMatrix("test_data/expression_matrix.txt")
  exprData <- ExpressionData(data)
  filteredList <- applyFilters(exprData, 0, TRUE, TRUE)
  filteredData <- filteredList[[1]]
  expect_equal(filteredData@data, data)
})

test_that("All filters reduce rows, not columns", {
  data <- readExprToMatrix("test_data/expression_matrix.txt")
  exprData <- ExpressionData(data)
  filteredList <- applyFilters(exprData, 20, FALSE, FALSE)
  filteredData <- filteredList[[1]]
  
  expect_equal(ncol(filteredData@thresholdFilter), ncol(data))
  expect_true(nrow(filteredData@thresholdFilter) <= nrow(data))
  
  expect_equal(ncol(filteredData@fanoFilter), ncol(data))
  expect_true(nrow(filteredData@fanoFilter) <= nrow(data))
  
})

