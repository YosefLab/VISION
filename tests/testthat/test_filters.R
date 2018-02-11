context("Filter Expression Data")

test_that("Threshold filter returns null on test matrix", {

  data <- matrix(c(c(1,0,0), c(0,0,1), c(0,0,0)), nrow=3, ncol=3, byrow=T)
  genes_passing <- filterGenesThreshold(data, 3)
  expect_equal(length(genes_passing), 0)
})

test_that("All filters reduce rows, not columns", {
  data <- readExprAsMatrix("test_data/expression_matrix.txt")
  exprData <- ExpressionData(data)
  filters <- c("fano")
  genes_passing <- applyFilters(getExprData(exprData), 20, filters)

  expect_true(length(genes_passing) <= nrow(data))
  expect_true(
      length(intersect(genes_passing, rownames(data))) ==
          length(genes_passing))
})
