setMethod("initialize", signature(.Object="Signature"),
          function(.Object, sigDict, name, source, metaData="") {
            if (missing(sigDict)) {
              stop("Missing sigDict information.")
            } else if (missing(name)) {
              stop("Missing signature name.")
            } else if (missing(source)) {
              stop("Missing source file.")
            }
            
            .Object@sigDict = sigDict
            .Object@name = name
            .Object@source = source
            .Object@metaData = metaData
            
            return(.Object);
            
          }
)

setMethod("sigEqual", signature(object="Signature"),
          function(object, compareSig) {
            if (class(compareSig) != "Signature") {
              stop("Input compareSig is not a Signature data type.")
            }
            
            if (object@name != compareSig@name) {
              return(FALSE)
            }
            
            return(object@sigDict == compareSig@sigDict)
            
          }
)

