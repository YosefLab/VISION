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
      data = "data.frame"
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
      debug = "logical",
      lean = "logical",
      subsample_size = "numeric",
      min_signature_genes = "numeric",
      qc = "logical",
      projections = "character",
      weights = "character",
      threshold = "numeric",
      sig_norm_method = "character",
      sig_score_method = "character",
      exprData = "data.frame",
      housekeepingData = "data.frame",
      sigData = "data.frame"),
    prototype = list(
      precomputed = NULL,
      output_dir = "FastProject_Output",
      nofilter = FALSE,
      nomodel = FALSE,
      pca_filter = FALSE,
      all_sigs = FALSE,
      debug = FALSE,
      lean = FALSE,
      subsample_size = 1000,
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

.Signature <- setClass("Signature",
    representation(
      sig_dict = "data.frame",
      signed = "logical",
      source = "character",
      name = "character"
))

.SignatureScores <- setClass("SignatureScores",
    representation(
      scores = "list",
      name = "character",
      sample_labels = "list",
      isFactor = "logical",
      isPrecomputed = "logical",
      numGenes = "numeric"
))
