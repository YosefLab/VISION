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
