###
# Wrapper classes for different types of data
#
# Wrapper classes for ExpressionData, ProbabilityData, and PCData
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

.ExpressionData <- setClass("ExpressionData", 
    representation(
      data = "data.frame",
      row_labels = "data.frame",
      col_labels = "data.frame"
))

.PCData <- setClass("PCData",
    representation(
      variance = "numeric",
      row_labels = "list",
      parent_data = "ExpressionData"
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
