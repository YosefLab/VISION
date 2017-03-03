require(cogena)
require(data.table)

readExprToMatrix <- function(filename, delimiter="\t") {
  
  message("Loading data from ", filename, " ...")
  
  data <- as.matrix(read.table(filename, sep= delimiter, header=TRUE, row.names=1))
  data2 <- unique(data)
  return(data2)
}

readHKGToMatrix <- function(filename, delimiter="\t") {
  message("Loading data from ", filename, "...")
  return(read.table(filename, sep=delimiter))
}


readSignaturesGmtToMatrix <- function(filename) {
  #TODO: ADD FUNCTIONALITY TO READ IN LIST OF FILENAMES FOR SIGNATURES
  
  message("Loading data from ", filename, " ...")
  inp <- gmt2list(filename)
  
  # Create a list to store all signatures
  sig_data = c()
  
  fp <- strsplit(filename, "/")
  f <- last(fp[[1]])
  
  # Create a dictionary to easily compute signature sign
  keys <- c("none" = 1.0,  'up' = 1.0, 'plus' = 1.0, 'down' = -1.0, 'dn' = -1.0, 'minus' = -1.0)
  
  for (sig in names(inp)) {
    genes <- c()
    values <- c()
    
    sig_parts = strsplit(sig, "_")
    key <- "none"
    
    if (is.element(tolower(last(sig_parts[[1]])), names(keys))) {
      key <- tolower(last(sig_parts[[1]]))
    }
    
    contains = FALSE
    
    #if (sig %in% sig_data[,1]) {
    #  contains = TRUE  
    #}
    
    for (k in inp[[sig]]) {
      elem <- strsplit(k, ",")
      if (length(elem[[1]]) > 1) {
        genes <- c(genes, elem[[1]][1]) 
        values <- c(values, as.numeric(elem[[1]][2])*keys[key])
      } else {
        genes <- c(genes, elem[[1]][1])
        values <- c(values, 1.0 * keys[key])
      }
    }
    
    if (contains) {
      sig_data[sig, 3] <- c(sig_data[sig,3], list(genes))
      sig_data[sig, 4] <- c(sig_data[sig,4], list(values))
    }
    
    else{

      names(values) <- genes
      newSig <- Signature(values, "None", f, sig)
      sig_data <- c(sig_data, newSig)

    }
  }
  
  return (sig_data)
}

