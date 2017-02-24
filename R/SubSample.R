#' This module handles the splitting and joining of the data if sub-sampling is enabled


splitSamples <- function(data, sampleSize) {
  # Set seed so outputs are reproduceable 
  set.seed(RANDOM_SEED)
  
  data <- sample(data)
  sub <- data[1:sampleSize,]
  holdout <- data[sampleSize:nrow(data),]
  #print(dim(holdout))
  
  return(list(holdout, sub))
  
  
  
}