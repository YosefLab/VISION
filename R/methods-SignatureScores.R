
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


#' Generate signature scores based on an input data.frame
#'
#' These are either factors or numeric values with known values for each sample
#' Can be conditions of interest, batch annotations, qc information, etc.
#'
#' @param df a data.frame object where each column is a cell signature. Names of
#' columns are signatue names. Rownames must correspond to appropriate samples
#'
#' @note Factor signatures and numeric signatures are treated differently. Make
#' sure that the class of each input column is the correct one.
SigScoresFromDataframe <- function(df) {
    sigScores <- c()
    for(i in 1:NCOL(df)) {
        sigScores <- c(sigScores, SignatureScores(df[,i],
                                              colnames(df)[i],
                                              rownames(df),
                                              is.factor(df[,i]),
                                              TRUE, 0))
    }
    return(sigScores)
}
