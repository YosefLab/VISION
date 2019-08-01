###
# Wrapper classes for different types of data
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

setClassUnion("numericORNULL", members = c("numeric", "NULL"))
setClassUnion("matrixORSparse",
    members = c("matrix", "dgeMatrix", "dgCMatrix", "dgTMatrix"))
setClassUnion("matrixORNULL", members = c("matrix", "NULL"))
setClassUnion("dataframeORNULL", members = c("data.frame", "NULL"))

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
    Consistency = "matrix",
    pValue = "matrix",
    FDR = "matrix",
    sigClusters = "list"
))

PCAnnotatorData <- setClass("PCAnnotatorData",
    slots = c(
        pearsonCorr = "matrix",
        pearsonCorrFeatures = "matrixORNULL"
    )
)

Trajectory <- setClass("Trajectory",
    slots = c(
        adjMat = "matrix", # MxM connectivity for milestones (w/ lengths)
        progressions = "data.frame" # position of cells between milestones
            # rownames: cell (character)
            # columns:  from (character), to (character), position (numeric, 0 to 1)
))

TrajectoryProjection <- setClass("TrajectoryProjection",
    slots = c(
        name = "character", # Name of projection
        vData = "matrix", # Mx2, Coordinates for milestones
        pData = "matrix", # Nx2, Coordiantes for cells
        adjMat = "matrix" # MxM, Connectivity for milestones
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
setClassUnion("PCAnnotatorDataOrNULL", members=c("PCAnnotatorData", "NULL"))

Vision <- setClass("Vision",
    slots = c(
        nomodel = "logical",
        projection_genes = "character",
        weights = "matrix",
        threshold = "numeric",
        sig_norm_method = "character",
        sig_score_method = "character",
        exprData = "matrixORSparse",
        featureBarcodeData = "matrixORSparse",
        unnormalizedData = "matrixORSparse",
        housekeepingData = "character",
        sigData = "list",
        metaData = "data.frame",
        perm_wPCA = "logical",
        pool = "logical",
        sigScores = "matrix",
        cellsPerPartition = "numeric",
        SigConsistencyScores = "ProjectionDataOrNULL",
        FeatureBarcodeConsistencyScores = "data.frame",
        ClusterFeatureBarcodeScores = "list",
        ClusterSigScores = "list",
        TrajectoryConsistencyScores = "ProjectionDataOrNULL",
        TrajectoryConsistencyScoresFeatures = "dataframeORNULL",
        PCAnnotatorData = "PCAnnotatorDataOrNULL",
        projection_methods = "character",
        Projections = "list",
        TrajectoryProjections = "list", # list of TrajectoryProjection
        SigGeneImportance = "list",
        pools = "list",
        inputProjections = "list",
        name = "character",
        num_neighbors = "numericORNULL",
        latentSpace = "matrix",
        latentTrajectory = "Trajectory",
        version = "numeric",
        selections = "list",
        params = "list"),
    prototype = list(
        nomodel = FALSE,
        weights = matrix(NA, 1, 1),
        threshold = 0,
        sig_norm_method = "znorm_rows",
        exprData = matrix(NA, 1, 1),
        featureBarcodeData = matrix(NA, 1, 1),
        unnormalizedData = matrix(NA, 1, 1),
        housekeepingData = character(),
        sigData = list(),
        metaData = data.frame(),
        perm_wPCA = FALSE,
        pool = FALSE,
        sigScores = matrix(NA, 1, 1),
        cellsPerPartition = 100,
        SigConsistencyScores = NULL,
        FeatureBarcodeConsistencyScores = data.frame(),
        ClusterSigScores = list(),
        ClusterFeatureBarcodeScores = list(),
        TrajectoryConsistencyScores = NULL,
        TrajectoryConsistencyScoresFeatures = NULL,
        PCAnnotatorData = NULL,
        projection_methods = character(),
        Projections = list(),
        TrajectoryProjections = list(),
        SigGeneImportance = list(),
        pools = list(),
        inputProjections = list(),
        name = "",
        num_neighbors = NULL, 
        latentSpace = matrix(NA, 1, 1),
        latentTrajectory = NULL,
        version = 1.11,
        selections = list(),
        params = list()
))
