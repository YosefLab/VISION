
applyFilters <- function(data, threshold, nofilter, lean) {
  
  if(nofilter) {
  
    data <- filterGenesNovar(data)

  } else {

    data <- filterGenesThreshold(data, threshold)
    
    if (!lean) {
      
      print("not lean!")
    }
    
  }
  
  col_labels <- data[1,]
  row_labels <- data[,1]
  
  return(list(data, row_labels, col_labels))
}

# Remove genes with 0 variance
filterGenesNovar <- function(data) {
  
  d <- data.frame(data)
  return (as.matrix(subset(d, apply(d[,-1], 1, var) != 0)))
  
}

filteGenesThreshold <- function(data, threshold) {
  
  d <- data.frame(data)
  return (as.matrix(subset(d, rowSums(d[,-1]) >= threshold)))
  
}