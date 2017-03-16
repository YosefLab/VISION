setMethod("initialize", signature(.Object="SignatureScores"),
          function(.Object, scores, name, sample_labels, 
                   isFactor, isPrecomputed, numGenes = "numeric") {
            
            .Object@scores = scores
            .Object@name = name
            .Object@sample_labels = sample_labels
            .Object@isFactor = isFactor
            .Object@isPrecomputed = isPrecomputed
            .Object@numGenes = numGenes 
            
            return(.Object)
          }
)