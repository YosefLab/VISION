#' Reads in an expression matrix from a file path to a matrix
#'	@param filename: path to the file to be read in
#'	@param delimiter: delimiter to use while reading in the file. Default is "\\t"
#'
#'  @returns: Matrix read in from filename.
#'  @examples
#'  expr <- readExprToMatrix("data/expression_matrix.txt")
readExprToMatrix <- function(filename, delimiter="\t") {

  message("Loading data from ", filename, " ...")

  #data <- as.matrix(read.table(filename, sep= delimiter, header=TRUE, row.names=1))
  data <- as.matrix(data.table::fread(filename, sep=delimiter, skip=1, drop=1))
  rnames <- data.table::fread(filename, sep=delimiter, skip=1, select=1)
  cnames <- read.table(filename, sep=delimiter, header=T, row.names=1,nrows=1)
  rownames(data) <- as.vector(t(rnames))
  colnames(data) <- colnames(cnames)

  # check if there are any repeats in samples or genes
  uniqrows <- unique(rownames(data))
  #uniqcols <- unique(colnames(data))

  data <- data[uniqrows,]
  #data <- data[, uniqcols]
  return(data)
}

#' Reads in a .txt file containing Housekeeping Gene names
#'
#' @param filename path to the file to be read in
#' @param delimiter delimiter to use while reading in the file. Default is tab.
#' @return Table containing the houseekeeping gene names
#' @examples
#' hkg <- readHKGToMatrix("data/housekeeping_genes.txt")
readHKGToMatrix <- function(filename, delimiter="\t") {

  message("Loading data from ", filename, "...")
  return(read.table(filename, sep=delimiter))
}

#' Reads in a list of signature input files.
#'
#' @importFrom cogena gmt2list
#' @param filenames a list of paths to signature input files
#' @return List of Signature Objects
#' @examples
#' sigs <- readSignaturesInput(c("sigs1.txt", "sigs2.gmt"))
readSignaturesInput <- function(filenames) {

  # Create a dictionary to easily compute signature sign
  keys <- c("none" = 1.0,  'up' = 1.0, 'plus' = 1.0, 'down' = -1.0, 'dn' = -1.0, 'minus' = -1.0, 'mius'=-1.0)


  # Create a list to store all signatures
  sig_data <- list()
  # Create list to store all relevant sig names
  sig_names <- c()

  for (filename in filenames) {
    message("Loading data from ", filename, " ...")
    fp <- unlist(strsplit(filename, "/"))
    f <- fp[[length(fp)]]
    fsplit <- unlist(strsplit(f, "[.]"))

    file_ext <- fsplit[[length(fsplit)]]

    if (file_ext == "gmt") {
      inp <- gmt2list(filename)

      for (sig in names(inp)) {
        genes <- c()
        values <- c()

        sig_parts = unlist(strsplit(sig, "_"))
        key <- "none"

        if (is.element(tolower(sig_parts[[length(sig_parts)]]), names(keys))) {
          key <- tolower(sig_parts[[length(sig_parts)]])
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
        newSig <- Signature(values, sig, f, "", cluster=0)
        sig_data <- c(sig_data, newSig)
        sig_names <- c(sig_names, sig)
      }

    } else if (file_ext == "txt") {
      sigList <- list()
      inp <- as.matrix(read.table(filename, sep="\t"))

      for (r in 1:nrow(inp)) {
        dat <- inp[r,]
        sigName <- dat[[1]]

        if (length(dat) == 2) {
          key <- "none"
          gene <- toupper(dat[[2]])
        } else {
          key <- tolower(dat[[2]])
          gene <- toupper(dat[[3]])
        }

        if (sigName %in% sig_names) {
          next
        }

        if (sigName %in% names(sigList)) {
          sigList[[sigName]][gene] <- as.numeric(1.0*keys[key])
        } else {
          sigList[[sigName]] <- c()
          sigList[[sigName]][[gene]] <- as.numeric(1.0*keys[key])
        }

      }

      for (sig in names(sigList)) {
        newSig <- Signature(sigList[[sig]], sig, f, "", isPrecomputed=FALSE, isFactor=FALSE, cluster=0)
        sig_data <- c(sig_data, newSig)
        sig_names <- c(sig_names, sig)
      }
    }
  }

  names(sig_data) <- sig_names

  return (sig_data)
}
#'  Reads precomputed signature values form a tab-delimited text file
#'
#'  First row of the file contains sample labels that the signatures correspond with
#'
#'  Each subsequent row contains a signature name in the first column,
#'  followed by the signature type (either 'numerical' or 'factor')
#'  followed by the signature values, one for each sample label in the file
#'
#'  @param  filename Filename to read the precomputed signatures from
#'  @param sampleLables List of labels for which we want the signature scores
#'  @param delimiter Delimiter to use when reading in file. Default is tab.
#'  @returns List of Signature data types
#'  @examples
#'  expr <- readExprToMatrix("expression.txt")
#'  samples <- colnames(expr)
#'  precomputed <- readPrecomputed("precomputed_sigs.txt", samples)
readPrecomputed <- function(filename, sampleLabels, delimiter="\t") {


  message("Loading data from ", filename)

  f <- as.matrix(read.table(filename, sep=delimiter))
  l1 <- f[1,]
  l2 <- f[2,]

  if ( (length(l1) != length(l2)) && (length(l1) != length(l2)-2)) {
    stop("Error in header line of precomputed signature file.\n
            First row should contain tab-separated list of samples")
  }

  # Make sure that l1 is in the same order as sample_labels; match indices between signatrues in file & sampleLabels
  sampleLabels <- as.character(as.matrix(sampleLabels))
  colnames(f) <- l1

  f[,3:ncol(f)] <- f[,sampleLabels]
  colnames(f) <- f[1,]

  sigScores <- c()

  for (s in 2:nrow(f)) {

    sigName <- f[s,1]
    sigType <- tolower(f[s,2])
    sigVals <- f[s,3:ncol(f)]
    if (sigType == "numerical") {
      sigIsFactor <- FALSE
      tryCatch(
        sigVals <- as.numeric(as.matrix(sigVals))
      , warning=function(w) {
        message("Error in precomputed signature:", sigName)
        ind <- which(is.na(sigVals))
        stop("Error in sample ", sample_labels[ind-2])
      })

    } else if (sigType == "factor") {
      sigIsFactor <- TRUE
      sigVals <- as.numeric(as.matrix(sigVals))

    } else {
      stop("Column 2 of precomputed signature file should specify either 'numerical' or 'factor'")
    }

    sigScores <- c(sigScores, SignatureScores(sigVals, sigName, sampleLabels, sigIsFactor, TRUE, 0))
  }

  return(sigScores)

}

