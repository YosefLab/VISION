context("Calculate Signature Background")

test_that("Can generate random background signatures", {

  data_file <- "test_data/expression_matrix.txt"
  sig_file <- "test_data/published_signatures/tcga_sigs.txt"

  object <- FastProject(data_file, sig_file)

  # Mock the weights
  ed <- getExprData(object@exprData)
  weights <- matrix(
                    rep(1, length(ed)),
                    nrow=nrow(ed)
                    )
  rownames(weights) <- rownames(ed)
  colnames(weights) <- colnames(ed)
  object@weights <- weights

  backgroundSigs <- calculateSignatureBackground(object, 10)

  expect_is(backgroundSigs, "list")
  expect_gt(length(backgroundSigs), 0)
  expect_equal(length(backgroundSigs) %% 10, 0)

  sample_sig <- backgroundSigs[[1]]

  expect_s4_class(sample_sig, "SignatureScores")

  expect_equal(length(sample_sig@scores), ncol(ed))
  expect_equal(length(sample_sig@sample_labels), ncol(ed))

  # Test this file too - has more signatures
  sig_file <- "test_data/published_signatures/h.all.v5.2.symbols.gmt"

  backgroundSigs <- calculateSignatureBackground(object, 10)

  expect_is(backgroundSigs, "list")
  expect_gt(length(backgroundSigs), 0)
  expect_equal(length(backgroundSigs) %% 10, 0)

  sample_sig <- backgroundSigs[[1]]

  expect_s4_class(sample_sig, "SignatureScores")

  expect_equal(length(sample_sig@scores), ncol(ed))
  expect_equal(length(sample_sig@sample_labels), ncol(ed))

})
