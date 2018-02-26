context("Calculate Signature Background")

test_that("Can generate random background signatures", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"

  object <- FastProject(data_file, sig_file)

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

  expect_is(backgroundSigs, "list")
  expect_gt(length(backgroundSigs), 0)

  expect_equal(length(backgroundSigs) %% 10, 0)

  sample_sig <- backgroundSigs[[1]]
  expect_is(sigAssignments, "factor")
  expect_equal(length(sigAssignments), length(object@sigData))

  expect_s4_class(sample_sig, "SignatureScores")

  expect_equal(length(sample_sig@scores), ncol(ed))
  expect_equal(names(sample_sig@scores), colnames(ed))

  # Test this file too - has more signatures
  sig_file <- "test_data/published_signatures/h.all.v5.2.symbols.gmt"

  out <- calculateSignatureBackground(object, 10)
  backgroundSigGroups <- out[[1]]
  sigAssignments <- out[[2]]

  expect_is(backgroundSigGroups, "list")
  expect_gt(length(backgroundSigGroups), 0)

  backgroundSigs <- backgroundSigGroups[[1]]

  expect_is(backgroundSigs, "list")
  expect_gt(length(backgroundSigs), 0)
  expect_equal(length(backgroundSigs) %% 10, 0)

  sample_sig <- backgroundSigs[[1]]

  expect_s4_class(sample_sig, "SignatureScores")

  expect_equal(length(sample_sig@scores), ncol(ed))
  expect_equal(names(sample_sig@scores), colnames(ed))

})
