###
# Wrapper classes for different types of data
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

setClassUnion('numericORNULL', members=c('numeric', 'NULL'))
setClassUnion('matrixORSparse', members=c("matrix", "dgCMatrix", "dgTMatrix"))

Cluster <- setClass("Cluster",
    slots = c(
    method = "character",
    param = "numeric",
    centers = "matrix",
    data = "matrixORSparse"
    )
)

ProjectionData <- setClass("ProjectionData",
    slots = c(
    sigProjMatrix = "matrix",
    pMatrix = "matrix",
    sigClusters = "list",
    emp_pMatrix = "matrix"
))

PCAnnotatorData <- setClass("PCAnnotatorData",
    slots = c(pearsonCorr = "matrix")
)

TreeProjection <- setClass("TreeProjection",
    slots = c(
    pData = "matrix",
    name = "character",
    vData = "matrix",
    adjMat = "matrix",
    edgeAssoc = "matrix",
    edgePos = "numeric"
))

TreeProjectionData <- setClass("TreeProjectionData",
    contains = c("ProjectionData"),
    slots = c(
        latentTree = "TreeProjection",
        projections = "list",
        treeScore = "numericORNULL"
))

ServerExpression <- setClass("ServerExpression",
    slots = c(
    data = "matrix",
    sample_labels = "character",
    gene_labels = "character"
))

ServerSigProjMatrix <- setClass("ServerSigProjMatrix",
    slots = c(
    zscores = "matrix",
    pvals = "matrix",
    proj_labels = "character",
    sig_labels = "character"
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
))

setClassUnion("ProjectionDataOrNULL", members=c("ProjectionData", "NULL"))
setClassUnion("TreeProjectionDataOrNULL", members=c("TreeProjectionData", "NULL"))
setClassUnion("PCAnnotatorDataOrNULL", members=c("PCAnnotatorData", "NULL"))

FastProject <- setClass("FastProject",
    slots = c(
        nomodel = "logical",
        projection_genes = "character",
        lean = "logical",
        weights = "matrix",
        threshold = "numeric",
        sig_norm_method = "character",
        sig_score_method = "character",
        trajectory_method = "character",
        exprData = "matrixORSparse",
        initialExprData = "matrixORSparse",
        unnormalizedData = "matrixORSparse",
        initialUnnormalizedData = "matrixORSparse",
        housekeepingData = "character",
        sigData = "list",
        metaData = "data.frame",
        initialMetaData = "data.frame",
        perm_wPCA = "logical",
        pool = "logical",
        sigScores = "matrix",
        cellsPerPartition = "numeric",
        SigConsistencyScores = "ProjectionDataOrNULL",
        ClusterSigScores = "list",
        TreeProjectionData = "TreeProjectionDataOrNULL",
        PCAnnotatorData = "PCAnnotatorDataOrNULL",
        Projections = "list",
        pools = "list",
        inputProjections = "list",
        name = "character",
        cluster_variable = "character",
        latentSpace = "matrix",
        initialLatentSpace = "matrix",
        version = "numeric"),
    prototype = list(
        nomodel = FALSE,
        lean = FALSE,
        weights = matrix(NA, 1, 1),
        threshold = 0,
        sig_norm_method = "znorm_rows",
        trajectory_method = "None",
        exprData = matrix(NA, 1, 1),
        initialExprData = matrix(NA, 1, 1),
        unnormalizedData = matrix(NA, 1, 1),
        initialUnnormalizedData = matrix(NA, 1, 1),
        housekeepingData = character(),
        sigData = list(),
        metaData = data.frame(),
        initialMetaData = data.frame(),
        perm_wPCA = FALSE,
        pool = FALSE,
        sigScores = matrix(NA, 1, 1),
        cellsPerPartition = 100,
        SigConsistencyScores = NULL,
        ClusterSigScores = list(),
        TreeProjectionData = NULL,
        PCAnnotatorData = NULL,
        Projections = list(),
        pools = list(),
        inputProjections = list(),
        name = "",
        cluster_variable = "",
        latentSpace = matrix(NA, 1, 1),
        initialLatentSpace = matrix(NA, 1, 1),
        version = 1.0
))
