context("Create Vision Object")

test_that("Can create the Vision Object", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"
  meta_file <- "test_data/precomputed_sigs.txt"
  meta <- read.table(meta_file, header = TRUE, row.names = 1)

  data <- read.table(data_file, sep = "\t", header = TRUE)

  expect_s4_class(
    Vision(data, sig_file),
    "Vision"
  )

  expect_s4_class(
    Vision(data, sig_file, meta = meta),
    "Vision"
  )

})
