require(cogena)

readTextToMatrix <- function(filename, delimiter="\t") {
  
  message("Loading data from ", filename, " ...")
  
  data <- data.frame(read.table(filename, sep= delimiter))
  
  data2 <- unique(data)

  return(data2)
}


readSignaturesGmtToMatrix <- function(filename) {
  #TODO: ADD FUNCTIONALITY TO READ IN LIST OF FILENAMES FOR SIGNATURES
  
  message("Loading data from ", filename, " ...")
  inp <- gmt2list(filename)
  
  header <- c("signature", "description", "genes", "expression_values", "file_of_origin")
  sig_data <- matrix(header, nrow=1, ncol=5)
  
  fp <- strsplit(filename, "/")
  f <- last(fp[[1]])
  
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
    
    if (sig %in% sig_data[,1]) {
      contains = TRUE  
    }
    
    for (k in inp[[sig]]) {
      elem <- strsplit(k, ",")
      if (length(elem[[1]]) > 1) {
        genes <- c(genes, elem[[1]][1])
        values <- c(values, as.numeric(elem[[1]][2]) * keys[key])
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
      entry <- c(sig,"None", list(genes), list(values), f)
      sig_data <- rbind(sig_data, matrix(entry, nrow=1, ncol=5))
      rn <- c(rownames(sig_data), sig)
    }
  }
  
  return (data.frame(sig_data))
}

