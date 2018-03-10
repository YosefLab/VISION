###
# Wrapper classes for different types of data
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

ProjectionData <- setClass("ProjectionData",
    slots = c(
    projections = "list",
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
    slots = c(pearsonCorr = "matrix")
)

Projection <- setClass("Projection",
    slots = c(
    name = "character",
    pData = "matrixORSparse",
    weights = "matrixORSparse"
))

TreeProjection <- setClass("TreeProjection",
    contains = "Projection",
    slots = c(
    vData = "matrix",
    adjMat = "matrix",
    edgeAssoc = "matrix",
    edgePos = "numeric"
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
    isMeta = "logical",
    isFactor = "logical"
    ),
    prototype = list(
    metaData = "",
    isMeta=FALSE,
    isFactor=FALSE
))

SignatureScores <- setClass("SignatureScores",
    slots = c(
    scores = "vector",
    name = "character",
    isFactor = "logical",
    isMeta = "logical",
    numGenes = "numeric"
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
        exprData = "matrixORSparse",
        initialExprData = "matrixORSparse",
        housekeepingData = "character",
        sigData = "list",
        metaData = "data.frame",
        initialMetaData = "data.frame",
        perm_wPCA = "logical",
        pool = "logical",
        sigScores = "list",
        cellsPerPartition = "numeric",
        ClusterProjectionData = "ProjectionData",
        ProjectionData = "ProjectionData",
        TreeProjectionData = "TreeProjectionData",
        PCAnnotatorData = "PCAnnotatorData",
        sigMatrix = "data.frame",
        pools = "list",
        inputProjections = "list",
        name = "character",
        cluster_variable = "character",
        latentSpace = "matrix",
        initialLatentSpace = "matrix"),
    prototype = list(
        nomodel = FALSE,
        pca_filter = FALSE,
        lean = FALSE,
        min_signature_genes = 5,
        weights = NULL,
        threshold = 0,
        sig_norm_method = "znorm_rows",
        trajectory_method = "None",
        exprData = NULL,
        housekeepingData = NULL,
        sigData = NULL,
        metaData = NULL,
        initialMetaData = NULL,
        perm_wPCA = FALSE,
        pool = FALSE,
        sigScores = NULL,
        cellsPerPartition = 100,
        ClusterProjectionData = NULL,
        ProjectionData = NULL,
        TreeProjectionData = NULL,
        PCAnnotatorData = NULL,
        sigMatrix = NULL,
        pools = NULL,
        inputProjections = NULL,
        name = "",
        cluster_variable = "",
        latentSpace = NULL,
        initialLatentSpace = NULL
))
