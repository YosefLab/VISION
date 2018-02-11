
#' Initializes a new SignatureScores object
#' Should not be called directly, instead use the `new` syntax
#'
#' @param scores Signature scores
#' @param name Name of signature
#' @param isFactor Indicates whether or not this is a factor signature
#' @param isMeta Indicates whether or not this score was derived from meta data
#' @param numGenes The number of genes used to calculate the score
#' @return New SignatureScores object
SignatureScores <- function(scores, name,
                    isFactor, isMeta, numGenes) {
    .Object <- new("SignatureScores", scores = scores, name = name,
                   isFactor = isFactor, isMeta = isMeta,
                   numGenes = numGenes)
    return(.Object)
    }


#' Generate signature scores based on an input data.frame
#'
#' These are either factors or numeric values with known values for each sample
#' Can be conditions of interest, batch annotations, qc information, etc.
#'
#' @param df a data.frame object where each column is a cell signature. Names of
#' columns are signatue names. Rownames must correspond to appropriate samples
#' @param sampleLabels a character vector representing sample names in
#' the expression matrix
#' @return a data.frame object where each column is a cell signature. Names of
#' columns are signatue names. Rownames must correspond to appropriate samples
#' @note Factor signatures and numeric signatures are treated differently. Make
#' sure that the class of each input column is the correct one.
SigScoresFromDataframe <- function(df, sampleLabels) {

    common <- intersect(row.names(df), sampleLabels)
    if (length(common) != length(sampleLabels)){
        stop("Provided meta data dataframe must have same sample labels as the expression matrix")
    }

    df <- df[sampleLabels, , drop = FALSE]

    return(df)
}
