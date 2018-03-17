#' @rdname FastProject-class
#' @export
setGeneric("FastProject", function(data, ...) {
    standardGeneric("FastProject")
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

setGeneric("computeKNNWeights", function(object, K) {
    standardGeneric("computeKNNWeights")
})
