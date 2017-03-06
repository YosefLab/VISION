###
# Wrapper classes for different types of data
#
# Wrapper classes for ExpressionData, ProbabilityData, and PCData
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

ExpressionData <- setClass("ExpressionData", 
    representation(
      data = "matrix"
))

FastProject <- setClass("FastProject",
    slots = c(
      data_file = "character",
      housekeeping = "character",
      signatures = "character",
      precomputed = "character",
      output_dir = "character",
      nofilter = "logical",
      nomodel = "logical",
      pca_filter = "logical",
      all_sigs = "logical",
      debug = "numeric",
      lean = "logical",
      subsample_size = "numeric",
      min_signature_genes = "numeric",
      qc = "logical",
      projections = "character",
      weights = "matrix",
      threshold = "numeric",
      sig_norm_method = "character",
      sig_score_method = "character",
      exprData = "matrix",
      housekeepingData = "data.frame",
      sigData = "list"),
    prototype = list(
      precomputed = NULL,
      output_dir = "FastProject_Output",
      nofilter = FALSE,
      nomodel = FALSE,
      pca_filter = FALSE,
      all_sigs = FALSE,
      debug = 0,
      lean = FALSE,
      subsample_size = 0,
      min_signature_genes = 5,
      qc=FALSE,
      projections = NULL,
      weights = NULL,
      threshold = 0,
      sig_norm_method = "znorm_rows",
      sig_score_method ="weighted avg",
      exprData = NULL,
      housekeepingData = NULL,
      sigData = NULL
))

Signature <- setClass("Signature",
    slots = c(
      sigDict = "vector",
      name = "character",
      source = "character",
      metaData = "character"
    ),
    prototype = list(
      metaData = ""
    )
)

SignatureScores <- setClass("SignatureScores",
    slots = c(
      scores = "vector",
      name = "character",
      sample_labels = "list",
      isFactor = "logical",
      isPrecomputed = "logical",
      numGenes = "numeric"
))
