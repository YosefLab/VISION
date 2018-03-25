context("Calculate Signature Background")

test_that("Can generate random background signatures", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"

  data <- read.table(data_file, sep = "\t", header = TRUE)

  object <- FastProject(data, sig_file)

  # Mock the weights
  ed <- object@exprData
  weights <- matrix(
                    rep(1, length(ed)),
                    nrow=nrow(ed)
                    )
  rownames(weights) <- rownames(ed)
  colnames(weights) <- colnames(ed)
  object@weights <- weights

  out <- calculateSignatureBackground(object, 10)
  backgroundSigGroups <- out[[1]]
  sigAssignments <- out[[2]]

  expect_is(backgroundSigGroups, "list")
  expect_gt(length(backgroundSigGroups), 0)

  backgroundSigs <- backgroundSigGroups[[1]]

  expect_is(backgroundSigs, "matrix")
  expect_gt(ncol(backgroundSigs), 0)

  expect_equal(ncol(backgroundSigs), 10)
  expect_equal(nrow(backgroundSigs), ncol(object@exprData))

  expect_is(sigAssignments, "factor")
  expect_equal(length(sigAssignments), length(object@sigData))

  # Test this file too - has more signatures
  sig_file <- "test_data/published_signatures/h.all.v5.2.symbols.gmt"

  out <- calculateSignatureBackground(object, 10)
  backgroundSigGroups <- out[[1]]
  sigAssignments <- out[[2]]

  expect_is(backgroundSigGroups, "list")
  expect_gt(length(backgroundSigGroups), 0)

  backgroundSigs <- backgroundSigGroups[[1]]

  expect_is(backgroundSigs, "matrix")
  expect_equal(ncol(backgroundSigs), 10)
  expect_equal(nrow(backgroundSigs), ncol(object@exprData))

})
