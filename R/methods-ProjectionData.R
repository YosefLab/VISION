#' Wrapper for storing all relevant information for a given consistency test
#'
#' Stores a list of consistency test results
#' Also stores clusters, a signature to projection matrix

#' Initializes a ProjectionData object for storing the consistency test results
#'
#' @param Consistency Matrix storing the consistency score per signature pair
#' @param pValue Matrix storing the p values for each projection, signature pair
#' @param FDR Matrix storing the FDR for each p-value (BH method)
#' @param sigClusters the signature clusters
#' @return ProjectionData object
ProjectionData <- function(Consistency, pValue, FDR, sigClusters) {

            .Object <- new("ProjectionData", Consistency = Consistency,
                           pValue = pValue, FDR = FDR, sigClusters = sigClusters)

            return(.Object)
            }
