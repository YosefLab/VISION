require(cogena)
require(reshape2)
require(data.table)

readTextToMatrix <- function(filename, delimiter) {
  
  sprintf("Loading data from %s ...", filename)
  
  data <- as.matrix(read.table(filename, sep= delimiter, header=TRUE, row.names=1))
  
  data2 <- unique(data)
  
  col_labels <- data2[1,]
  row_labels <- data2[,1]
  
  #if (nrows(data2) != nrows(data)) {
  #  printf("WARNING: Row labels (gene identifiers) are not unique.\n")
    #, paste(row_labels_copy))
  #  printf("Removing these genes...")
  #}
  
  return(list(data2, col_labels, row_labels))
}


## TODO: find a file with a "_down" or "_up" and implement this featurization ##
readSignaturesGmtToMatrix <- function(filename) {

  inp <- gmt2list(filename)
  
  header <- c("signature", "description", "genes", "expression_values", "file_of_origin")
  sig_data <- matrix(header, nrow=1, ncol=5)
  
  fp <- strsplit(filename, "/")
  f <- last(fp[[1]])
  
  for (sig in names(inp)) {
    genes <- c()
    values <- c()
    for (k in inp[[sig]]) {
      elem <- strsplit(k, ",")
      if (length(elem[[1]]) > 1) {
        genes <- c(genes, elem[[1]][1])
        values <- c(values, as.numeric(elem[[1]][2]))
      } else {
        genes <- c(genes, elem[[1]][1])
        values <- c(values, 1.0)
      }
    }
    entry <- c(sig, "None", list(genes), list(values), f)
    sig_data <- rbind(sig_data, matrix(entry, nrow=1, ncol=5))
  }
  
  return(as.matrix(sig_data, row.names=1, header=TRUE))

}

