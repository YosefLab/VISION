###
# Wrapper classes for different types of data
#
# This was created to organize nuances in how signature
# scores, distance matrices, and anything else, is computed
# on the different types of data.

setClassUnion("numericORNULL", members = c("numeric", "NULL"))
setClassUnion("matrixORSparse", members = c("matrix", "dgCMatrix"))
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

NormData <- setClass("NormData",
    slots = c(
    colOffsets = "numeric",
    colScaleFactors = "numeric",
    rowOffsets = "numeric",
    rowScaleFactors = "numeric",
    data = "matrixORSparse"
    ),
    validity = function(obj){
        isValid <- nrow(obj@data) == length(obj@rowOffsets)
        isValid <- isValid && (nrow(obj@data) == length(obj@rowScaleFactors))
        isValid <- isValid && (ncol(obj@data) == length(obj@colOffsets))
        isValid <- isValid && (ncol(obj@data) == length(obj@colScaleFactors))
        return(isValid)
    },
)

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

setClassUnion("PCAnnotatorDataOrNULL", members=c("PCAnnotatorData", "NULL"))

Vision <- setClass("Vision",
    slots = c(
        exprData = "matrixORSparse",
        featureBarcodeData = "matrixORSparse",
        unnormalizedData = "matrixORSparse",
        sigData = "list",
        metaData = "data.frame",
        SigScores = "matrix",
        LocalAutocorrelation = "list",
        TrajectoryAutocorrelation = "list",
        ClusterComparisons = "list",
        PCAnnotatorData = "PCAnnotatorDataOrNULL",
        Projections = "list",
        TrajectoryProjections = "list", # list of TrajectoryProjection
        SigGeneImportance = "list",
        Pools = "list",
        LatentSpace = "matrix",
        LatentTrajectory = "Trajectory",
        Viewer = "list",
        params = "list",
        version = "numeric"
        ),
    prototype = list(
        exprData = matrix(NA, 1, 1),
        featureBarcodeData = matrix(NA, 1, 1),
        unnormalizedData = matrix(NA, 1, 1),
        sigData = list(),
        metaData = data.frame(),
        SigScores = matrix(NA, 1, 1),
        LocalAutocorrelation = list(),
        TrajectoryAutocorrelation = list(),
        ClusterComparisons = list(),
        PCAnnotatorData = NULL,
        Projections = list(),
        TrajectoryProjections = list(),
        SigGeneImportance = list(),
        Pools = list(),
        LatentSpace = matrix(NA, 1, 1),
        LatentTrajectory = NULL,
        Viewer = list(),
        params = list(),
        version = 1.11
))
