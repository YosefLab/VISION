context("Load Precomputed Signatures")

test_that("Can read precomputed signatures from a file.", {

  data_file <- "test_data/expression_matrix.txt"
  sample_labels <- colnames(
                            read.table(data_file)
                            )

  precomp_file <- "test_data/precomputed_sigs.txt"

  sigs <- readPrecomputed(precomp_file, sample_labels)

  expect_equal(length(sigs), 3)
})

test_that("Can load precomputed signatures from a dataframe", {
    sample_labels <- c("A", "B", "C", "D")
    time <- c(1, 2, 1, 2)
    batch <- as.factor(c(1, 1, 2, 2))
    donor <- as.factor(c("A", "A", "B", "B"))

    df <- data.frame(
                    time=time,
                    batch=batch,
                    donor=donor,
                    row.names=sample_labels
                    )

    sigs <- SigScoresFromDataframe(df, sample_labels)
    expect_equal(length(sigs), 3)

    # provided precomputed df can be a super-set of the samples in the matrix
    sample_labels_ok = c("A", "B", "C")
    sigs <- SigScoresFromDataframe(df, sample_labels_ok)
    expect_equal(length(sigs), 3)

    # Must have all the samples that are in the df
    sample_labels_bad = c("A", "B", "C", "D", "E")
    expect_error(SigScoresFromDataframe(df, sample_labels_bad))
  
})

