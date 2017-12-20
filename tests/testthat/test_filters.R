context("Filter Expression Data")

test_that("Threshold filter returns null on test matrix", {

  data <- matrix(c(c(1,0,0), c(0,0,1), c(0,0,0)), nrow=3, ncol=3, byrow=T)
  tdata <- filterGenesThreshold(data, 3)
  expect_equal(nrow(tdata), 0)
  expect_equal(ncol(tdata), 3)
})

test_that("All filters reduce rows, not columns", {
  data <- readExprAsMatrix("test_data/expression_matrix.txt")
  exprData <- ExpressionData(data)
  filters <- c("fano")
  filteredData <- applyFilters(exprData, 20, filters)
  
  expect_equal(ncol(filteredData@fanoFilter), ncol(data))
  expect_true(nrow(filteredData@fanoFilter) <= nrow(data))
  
})

