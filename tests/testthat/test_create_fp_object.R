context("Create FastProject Object")

test_that("Can create the FastProject Object", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"
  precomp_file <- "test_data/precomputed_sigs.txt"

  expect_s4_class(
    FastProject(data_file, sig_file),
    "FastProject"
  )

  expect_s4_class(
    FastProject(data_file, sig_file, precomputed=precomp_file),
    "FastProject"
  )

})
