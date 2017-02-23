setGeneric("Analyze", function(object) {
  standardGeneric("Analyze")
})

setGeneric("createOutputDirectory", function(object) {
  standardGeneric("createOutputDirectory")
})

setGeneric("getExprData", function(object) {
  standardGeneric("getExprData")
})

setGeneric("readExprData", function(object) {
  standardGeneric("readExprData")
})

setGeneric("updateExprData", function(object, newData=data.frame()) {
  standardGeneric("updateExprData")
})