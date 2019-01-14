#' @rdname VISION-class
#' @export
setGeneric("Vision", function(data, ...) {
    standardGeneric("Vision")
})

setGeneric("analyze", function(object, ...) {
    standardGeneric("analyze")
})

setGeneric("addProjection", function(object, ...) {
    standardGeneric("addProjection")
})

setGeneric("cluster", function(object, method, param) {
    standardGeneric("cluster")
})

setGeneric("saveAndViewResults", function(object, ...) {
    standardGeneric("saveAndViewResults")
})

setGeneric("updateProjection", function(object, name, data) {
    standardGeneric("updateProjection")
})

setGeneric("viewResults", function(object, ...) {
    standardGeneric("viewResults")
})

setGeneric("getSelections", function(object, ...) {
    standardGeneric("getSelections")
})

setGeneric("computeKNNWeights", function(object, ...) {
    standardGeneric("computeKNNWeights")
})

setGeneric("getProjections", function(object, ...) {
    standardGeneric("getProjections")
})

setGeneric("getLatentSpace", function(object, ...) {
    standardGeneric("getLatentSpace")
})

setGeneric("getLatentTrajectory", function(object, ...) {
    standardGeneric("getLatentTrajectory")
})

setGeneric("getSignatureScores", function(object, ...) {
    standardGeneric("getSignatureScores")
})

setGeneric("getSignatureConsistency", function(object, ...) {
    standardGeneric("getSignatureConsistency")
})

setGeneric("getSignatureDifferential", function(object, ...) {
    standardGeneric("getSignatureDifferential")
})

setGeneric("getMetaConsistency", function(object, ...) {
    standardGeneric("getMetaConsistency")
})

setGeneric("getMetaDifferential", function(object, ...) {
    standardGeneric("getMetaDifferential")
})
