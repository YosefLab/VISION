setGeneric("addSigData", function(object, data) {
  standardGeneric("addSigData")
})


setGeneric("Analyze", function(object) {
  standardGeneric("Analyze")
})

setGeneric("createOutputDirectory", function(object) {
  standardGeneric("createOutputDirectory")
})

setGeneric("generateOutputReport", function(fpout, url="http://127.0.0.1:8080/html/Results.html") {
  standardGeneric("generateOutputReport")
})

setGeneric("getExprData", function(object) {
  standardGeneric("getExprData")
})

setGeneric("getNormalizedCopy", function(object, func) {
  standardGeneric("getNormalizedCopy")
})

setGeneric("sigEqual", function(object, compareSig) {
  standardGeneric("sigEqual")
})

setGeneric("readExprData", function(object) {
  standardGeneric("readExprData")
})

setGeneric("updateExprData", function(object, newData) {
  standardGeneric("updateExprData")
})


setGeneric("updateProjection", function(object, name, data) {
  standardGeneric("updateProjection")
})

setGeneric("cluster", function(object, method, param) {
  standardGeneric("cluster")
})
