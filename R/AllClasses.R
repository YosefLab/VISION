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

# setClassUnion("treeorNull", members=c("phylo", "NULL"))
# setClassUnion("pythonorNull", members = c("python.builtin.object", "NULL"))

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
    data = "Matrix"
    ),
    validity = function(object){
        isValid <- nrow(object@data) == length(object@rowOffsets)

        isValid <- isValid && (
            nrow(object@data) == length(object@rowScaleFactors))

        isValid <- isValid && (
            ncol(object@data) == length(object@colOffsets))

        isValid <- isValid && (
            ncol(object@data) == length(object@colScaleFactors))

        return(isValid)
    },
)

LCAnnotatorData <- setClass("LCAnnotatorData",
    slots = c(
        pearsonCorr = "matrix",
        pearsonCorrProteins = "matrixORNULL"
    )
)

Trajectory <- setClass("Trajectory",
    slots = c(
        adjMat = "matrix", # MxM connectivity for milestones (w/ lengths)
        progressions = "data.frame" # position of cells between milestones
            # rownames: cell (character)
            # columns:  from (character), to (character), position (numeric, 0 to 1)
))

setClassUnion("trajectoryORNULL", members = c("Trajectory", "NULL"))

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

setClassUnion("LCAnnotatorDataOrNULL", members=c("LCAnnotatorData", "NULL"))

Vision <- setClass("Vision",
    slots = c(
        exprData = "matrixORSparse",
        proteinData = "matrixORSparse",
        unnormalizedData = "matrixORSparse",
        sigData = "list",
        metaData = "data.frame",
        modData = "list",
        SigScores = "matrix",
        ModScores = "matrix",
        LocalAutocorrelation = "list",
        TrajectoryAutocorrelation = "list",
        ClusterComparisons = "list",
        LCAnnotatorData = "LCAnnotatorDataOrNULL",
        Projections = "list",
        TrajectoryProjections = "list", # list of TrajectoryProjection
        SigGeneImportance = "list",
        ModGeneImportance = "list",
        Pools = "list",
        LatentSpace = "matrix",
        LatentTrajectory = "trajectoryORNULL",
        Hotspot = "list",
        ModuleSignatureEnrichment = "list",
        ModuleHotspotScores = "data.frame",
        Viewer = "list",
        params = "list",
        version = "numeric"
        ),
    prototype = list(
        exprData = matrix(NA, 1, 1),
        proteinData = matrix(NA, 1, 1),
        unnormalizedData = matrix(NA, 1, 1),
        sigData = list(),
        metaData = data.frame(),
        modData=list(),
        SigScores = matrix(NA, 1, 1),
        ModScores = matrix(NA, 1, 1),
        LocalAutocorrelation = list(),
        TrajectoryAutocorrelation = list(),
        ClusterComparisons = list(),
        LCAnnotatorData = NULL,
        Projections = list(),
        TrajectoryProjections = list(),
        SigGeneImportance = list(),
        ModGeneImportance = list(),
        Pools = list(),
        LatentSpace = matrix(NA, 1, 1),
        LatentTrajectory = NULL,
        Hotspot = list(),
        ModuleSignatureEnrichment = list(),
        ModuleHotspotScores = data.frame(),
        Viewer = list(),
        params = list(),
        version = 1.2
))

# This is a hack to try and get R to play nice with the phylo class
# Which the Ape package doesn't expose
phylo <- setClass("phylo")

PhyloVision <- setClass("PhyloVision", contains = "Vision", 
  slots = c(tree = "phylo"),
  prototype = list(tree=NULL)
)
