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


readSignaturesInput <- function(filenames) {

  # Create a dictionary to easily compute signature sign
  keys <- c("none" = 1.0,  'up' = 1.0, 'plus' = 1.0, 'down' = -1.0, 'dn' = -1.0, 'minus' = -1.0)
  
  # Create a list to store all signatures
  sig_data <- c()
  # Create list to store all relevant sig names
  sig_names <- c()
  
  for (filename in filenames) {
    message("Loading data from ", filename, " ...")
    fp <- strsplit(filename, "/")
    f <- last(fp[[1]])
    
    if (file_ext(f) == "gmt") {
      inp <- gmt2list(filename)
      
      for (sig in names(inp)) {
        genes <- c()
        values <- c()
        
        sig_parts = strsplit(sig, "_")
        key <- "none"
        
        if (is.element(tolower(last(sig_parts[[1]])), names(keys))) {
          key <- tolower(last(sig_parts[[1]]))
        }
        
        if (sig %in% sig_names) {
          next; 
        }
        
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
    
        names(values) <- genes
        newSig <- Signature(values, sig, f, "")
        sig_data <- c(sig_data, newSig)
        sig_names <- c(sig_names, sig)
      }
      
    } else if (file_ext(filename) == "txt") {
      sigList <- list()
      inp <- as.matrix(read.table(filename, sep="\t"))
      
      for (r in 1:nrow(inp)) {
        dat <- inp[r,]
        sigName <- dat[[1]]

        if (length(dat) == 2) {
          key <- "none"
          gene <- dat[[2]]
        } else {
          key <- tolower(dat[[2]])
          gene <- dat[[3]]
        }
        
        if (sigName %in% sig_names) {
          next
        }
        
        if (sigName %in% names(sigList)) {
          sigList[[sigName]][gene] <- 1.0*keys[key]
        } else {
          sigList[[sigName]] <- list()
          sigList[[sigName]][[gene]] <- 1.0*keys[key]
        }
      
      }
      
      for (sig in names(sigList)) {
        newSig <- Signature(sigList[[sig]], sig, f, "")
        sig_data <- c(sig_data, newSig)
        sig_names <- c(sig_names, sig)
      } 
    }
  }
  
  return (sig_data)
}

