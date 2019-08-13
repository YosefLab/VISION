context("Calculate Signature Background")

test_that("Can generate random background signatures", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"

  data <- read.table(data_file, sep = "\t", header = TRUE)

  object <- Vision(data, sig_file)

  # Mock the weights
  ed <- object@exprData

  out <- generatePermutationNull(object@exprData, object@sigData, 10)
  randomSigs <- out$randomSigs
  sigAssignments <- out$sigAssignments
  randomSigAssignments <- out$randomSigAssignments

  expect_is(randomSigs, "list")
  expect_gt(length(randomSigs), 0)

  expect_is(sigAssignments, "factor")
  expect_equal(length(sigAssignments), length(object@sigData))
  expect_equal(
      length(levels(sigAssignments)),
      length(randomSigs) / 10
  )

  expect_is(randomSigAssignments, "factor")
  expect_equal(length(randomSigAssignments), length(randomSigs))
  expect_equal(
      length(levels(randomSigAssignments)),
      length(randomSigs) / 10
  )


  # Test this file too - has more signatures
  sig_file <- "test_data/published_signatures/h.all.v5.2.symbols.gmt"

  out <- generatePermutationNull(object@exprData, object@sigData, 10)
  randomSigs <- out$randomSigs
  sigAssignments <- out$sigAssignments
  randomSigAssignments <- out$randomSigAssignments

  expect_is(randomSigs, "list")
  expect_gt(length(randomSigs), 0)

  expect_is(sigAssignments, "factor")
  expect_equal(length(sigAssignments), length(object@sigData))
  expect_equal(
      length(levels(sigAssignments)),
      length(randomSigs) / 10
  )

  expect_is(randomSigAssignments, "factor")
  expect_equal(length(randomSigAssignments), length(randomSigs))
  expect_equal(
      length(levels(randomSigAssignments)),
      length(randomSigs) / 10
  )

})
