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

setGeneric("saveAndViewResults", function(fpout, ...) {
    standardGeneric("saveAndViewResults")
})

setGeneric("updateProjection", function(object, name, data) {
    standardGeneric("updateProjection")
})

setGeneric("viewResults", function(object, ...) {
    standardGeneric("viewResults")
})

setGeneric("computeKNNWeights", function(object, ...) {
    standardGeneric("computeKNNWeights")
})
