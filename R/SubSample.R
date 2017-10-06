#' This module handles the splitting and joining of the data if sub-sampling is enabled

#' randomally subsample the data
#' @param data the data to subsample
#' @param sampleSize the size of the subsample
#' @return a list:
#' \itemize{
#'     \item the data not included in the subsample
#'     \item the data included in the subsample
#' }
splitSamples <- function(data, sampleSize) {
    # Set seed so outputs are reproduceable
    set.seed(RANDOM_SEED)

    data <- sample(data)
    sub <- data[1:sampleSize,]
    holdout <- data[sampleSize:nrow(data),]

    return(list(holdout, sub))



}
