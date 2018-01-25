###
# Wrapper classes for different types of data
#
# Wrapper classes for ExpressionData, ProbabilityData, and PCData
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
setClassUnion('matrixORSparse', members=c("matrix", "data.frame", "dgCMatrix", "dgTMatrix"))

Cluster <- setClass("Cluster",
    slots = c(
    method = "character",
    param = "numeric",
    centers = "matrix",
    data = "matrixORSparse"
    )
)

ExpressionData <- setClass("ExpressionData",
    slots = c(
        data = "matrixORSparse",
        fanoFilter = "matrixORSparse"
))

FastProject <- setClass("FastProject",
    slots = c(
        nomodel = "logical",
        projection_genes = "character",
        lean = "logical",
        min_signature_genes = "numeric",
        weights = "matrix",
        threshold = "numeric",
        sig_norm_method = "character",
        sig_score_method = "character",
        trajectory_method = "character",
        exprData = "ExpressionData",
        allData = "matrixORSparse",
        housekeepingData = "character",
        sigData = "list",
        precomputedData= "list",
        perm_wPCA = "logical",
        pool = "logical",
        sigScores = "list",
        cellsPerPartition= "numeric",
        filterModuleList = "list",
        sigMatrix = "data.frame",
        pools = "list",
        inputProjections = "list",
        name = "character"),
    prototype = list(
        nomodel = FALSE,
        pca_filter = FALSE,
        lean = FALSE,
        min_signature_genes = 5,
        weights = NULL,
        threshold = 0,
        sig_norm_method = "znorm_rows",
        trajectory_method ="None",
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
        pools=NULL,
        inputProjections=NULL,
        name = ""
))

ProjectionData <- setClass("ProjectionData",
    slots = c(
    projections = "list",
    keys = "character", # names of projections
    sigProjMatrix = "matrix",
    pMatrix="matrix",
    sigClusters = "list",
    emp_pMatrix = "matrix"
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

setClassUnion('TreeProjDataORNULL', members=c('TreeProjectionData', 'NULL'))
FilterModuleData <- setClass("FilterModuleData",
    slots = c(
    filter = "character",
    genes = "character",
    ProjectionData = "ProjectionData",
    TreeProjectionData = "TreeProjDataORNULL",
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
    isFactor = "logical",
    isPrecomputed = "logical",
    numGenes = "numeric"
))
