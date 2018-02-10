context("Create FastProject Object")

test_that("Can create the FastProject Object", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"
  meta_file <- "test_data/precomputed_sigs.txt"
  meta <- read.table(meta_file, header = TRUE, row.names = 1)

  expect_s4_class(
    FastProject(data_file, sig_file),
    "FastProject"
  )

  expect_s4_class(
    FastProject(data_file, sig_file, meta = meta),
    "FastProject"
  )

})
