
#' Initializes a new SignatureScores object
#' Should not be called directly, instead use the `new` syntax
#'
#' @param .Object an object
#' @param scores Signature scores
#' @param name Name of signature
#' @param sample_labels Sample names of expression matrix
#' @param isFactor Indicates whether or not this is a factor signature
#' @param isPrecomputed Indicates whether or not this score was precomputed
#' @param numGenes The number of genes used to calculate the score
#' @return New SignatureScores object
setMethod("initialize", signature(.Object="SignatureScores"),
            function(.Object, scores, name, sample_labels,
                    isFactor, isPrecomputed, numGenes) {

            .Object@scores = scores
            .Object@name = name
            .Object@sample_labels = sample_labels
            .Object@isFactor = isFactor
            .Object@isPrecomputed = isPrecomputed
            .Object@numGenes = numGenes

            return(.Object)
            }
)
