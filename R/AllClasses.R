###
# Wrapper classes for different types of data
#
# Wrapper classes for ExpressionData, ProbabilityData, and PCData
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

setClassUnion('numericORNULL', members=c('numeric', 'NULL'))

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
        nofilter = "logical",
        nomodel = "logical",
        filters = "character",
        lean = "logical",
        min_signature_genes = "numeric",
        qc = "logical",
        projections = "character",
        weights = "matrix",
        threshold = "numeric",
        sig_norm_method = "character",
        sig_score_method = "character",
        exprData = "ExpressionData",
        allData = "matrix",
        housekeepingData = "character",
        sigData = "list",
        precomputedData= "list",
        perm_wPCA = "logical",
        pool = "logical",
        sigScores = "list",
        cellsPerPartition= "numeric",
        filterModuleList = "list",
        sigMatrix = "data.frame",
        fpParams = "list",
        pools = "list"),
    prototype = list(
        nofilter = FALSE,
        nomodel = FALSE,
        pca_filter = FALSE,
        lean = FALSE,
        min_signature_genes = 5,
        qc = FALSE,
        projections = NULL,
        weights = NULL,
        threshold = 0,
        sig_norm_method = "znorm_rows",
        sig_score_method ="weighted avg",
        exprData = NULL,
        housekeepingData = NULL,
        sigData = NULL,
        precomputedData = NULL,
        perm_wPCA = FALSE,
        pool = FALSE,
        sigScores = NULL,
        cellsPerPartition = 100,
        filterModuleList = NULL,
        sigMatrix = NULL,
        fpParams = NULL,
        pools=NULL
))

ProjectionData <- setClass("ProjectionData",
    slots = c(
    projections = "list",
    keys = "character", # names of projections
    sigProjMatrix = "matrix",
    pMatrix="matrix",
    sigClusters = "list"
))

TreeProjectionData <- setClass("TreeProjectionData",
    contains = c("ProjectionData"),
    slots = c(
    treeScore = "numericORNULL"
    # adjMat = "matrix"
))

PCAnnotatorData <- setClass("PCAnnotatorData",
    slots = c(
    fullPCA = "matrix",
    pearsonCorr = "matrix",
    loadings = "matrix"
))

Projection <- setClass("Projection",
    slots = c(
    name = "character",
    pData = "matrix",
    weights = "matrix"
))

TreeProjection <- setClass("TreeProjection",
    contains = "Projection",
    slots = c(
    vData = "matrix",
    adjMat = "matrix",
    edgeAssoc = "matrix",
    edgePos = "numeric"
))

FilterModuleData <- setClass("FilterModuleData",
    slots = c(
    filter = "character",
    genes = "character",
    ProjectionData = "ProjectionData",
    TreeProjectionData = "TreeProjectionData",
    PCAnnotatorData = "PCAnnotatorData"
))

ServerExpression <- setClass("ServerExpression",
    slots = c(
    data = "matrix",
    sample_labels = "character",
    gene_labels = "character"
))

ServerPCorr <- setClass("ServerPCorr",
    slots = c(
    data = "matrix",
    proj_labels = "character",
    sig_labels = "character"
))

ServerSigProjMatrix <- setClass("ServerSigProjMatrix",
    slots = c(
    data = "matrix",
    proj_labels = "character",
    sig_labels = "character"
))

ServerPMatrix <- setClass("ServerPMatrix",
    slots = c(
    data = "matrix",
    proj_labels = "character",
    sig_labels = "character"
))

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
))

SignatureScores <- setClass("SignatureScores",
    slots = c(
    scores = "vector",
    name = "character",
    sample_labels = "character",
    isFactor = "logical",
    isPrecomputed = "logical",
    numGenes = "numeric"
))
