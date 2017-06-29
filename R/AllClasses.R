###
# Wrapper classes for different types of data
#
# Wrapper classes for ExpressionData, ProbabilityData, and PCData
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

Cluster <- setClass("Cluster",
  slots = c(
    method = "character",
    param = "numeric",
    centers = "matrix",
    data = "matrix"
  )
)

ExpressionData <- setClass("ExpressionData", 
    slots = c(
      data = "matrix",
      fanoFilter = "matrix",
      noVarFilter = "matrix",
      thresholdFilter = "matrix"
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
      filters = "character",
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
      sigData = "list",
      precomputedData= "list",
      perm_wPCA = "logical",
      numCores = "numeric",
      approximate = "logical",
      optClust = "numeric"), 
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
      sigData = NULL,
      precomputedData=NULL,
      perm_wPCA=FALSE,
      numCores = 0,
      optClust = 0,
      approximate=F
))

FastProjectOutput <- setClass("FastProjectOutput",
    slots = c(
      exprData = "ExpressionData",
      projData = "list", 
      sigMatrix = "matrix",
      sigList = "list",
      sigClusters = "list"
    )
)

Projection <- setClass("Projection", 
    slots = c(
      name = "character",
      pData = "matrix"
    )
)

ProjectionData <- setClass("ProjectionData", 
    slots = c(
      filter = 'character', 
      projections = "list",
      genes = "character",
      keys = "character",
      sigProjMatrix = "matrix",
      pMatrix="matrix",
      PPT = "list"
    )
)

ServerExpression <- setClass("ServerExpression", 
    slots = c(
      data = "matrix", 
      sample_labels = "character",
      gene_labels = "character"
    )
)

ServerSigProjMatrix <- setClass("ServerSigProjMatrix",
    slots = c(
      data = "matrix",
      proj_labels = "character",
      sig_labels = "character"
    )
)

ServerPMatrix <- setClass("ServerPMatrix",
    slots = c(
      data = "matrix",
      proj_labels = "character",
      sig_labels = "character"
    )
)

Signature <- setClass("Signature",
    slots = c(
      sigDict = "vector",
      name = "character",
      source = "character",
      metaData = "character",
      isPrecomputed = "logical",
      isFactor = "logical",
      cluster = "numeric"
    ),
    prototype = list(
      metaData = "",
      isPrecomputed=FALSE,
      isFactor=FALSE,
      cluster=0
    )
)

SignatureScores <- setClass("SignatureScores",
    slots = c(
      scores = "vector",
      name = "character",
      sample_labels = "character",
      isFactor = "logical",
      isPrecomputed = "logical",
      numGenes = "numeric"
))
