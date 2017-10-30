#' Reads in an expression matrix from a file path to a matrix
#' @param filename path to the file to be read in
#' @param sep seperator to use while reading in the file. Default is "\\t"
#'
#' @return Matrix read in from filename.
readExprAsMatrix <- function(filename, sep="\t") {

    message("Loading data from ", filename, " ...")
    f <- basename(filename)
    fsplit <- unlist(strsplit(f, "[.]"))
    if (fsplit[[length(fsplit)]] == "rds") {
        return(readRDS(filename))
    }
    data <- as.matrix(read.table(filename, header = TRUE, sep=sep,
                                 row.names = 1))

    # make sure gene names are all upper case to facilitate easy matching later
    rownames(data) <- toupper(rownames(data))

    # check if there are any repeats in samples or genes
    data <- data[unique(rownames(data)),]

    return(data)
}

#' Reads in a list of signature input files.
#' @importFrom utils read.table
#' @importFrom cogena gmt2list
#' @param filenames a list of paths to signature input files
#' @return List of Signature Objects
readSignaturesInput <- function(filenames) {

    # Create a dictionary to easily compute signature sign
    keys <- c("none" = 1.0,  'up' = 1.0, 'plus' = 1.0, 'down' = -1.0,
            'dn' = -1.0, 'minus' = -1.0, 'mius'=-1.0)


    # Create a list to store all signatures
    sig_data <- list()

    for (filename in filenames) {
    message("Loading data from ", filename, " ...")

    fsplit <- strsplit(basename(filename), "\\.")[[1]]
    file_ext <- fsplit[2]

    if (file_ext == "gmt") {
        inp <- gmt2list(filename)

        for (sig in names(inp)) {
        genes <- c()
        values <- c()

        sig_parts = unlist(strsplit(sig, "_"))
        key <- "none"

        if (is.element(tolower(sig_parts[[length(sig_parts)]]), names(keys))) {
            key <- tolower(sig_parts[[length(sig_parts)]])
            sig_name <- paste(sig_parts[1:length(sig_parts)-1], collapse="_")
        } else {
            sig_name <- sig
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
        if (sig_name %in% names(sig_data)) {
            existingSig <- sig_data[[sig_name]]
            existingSig@sigDict[names(values)] <- values
            sig_data[[sig_name]] <- existingSig
        } else {
            newSig <- Signature(values, sig_name, filename, "", cluster=0)
            sig_data[[sig_name]] <- newSig
        }
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

        if (sigName %in% names(sig_data)) {
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
        newSig <- Signature(sigList[[sig]], sig, filename, "", isPrecomputed=FALSE,
                            isFactor=FALSE, cluster=0)
        sig_data[[sig]] <- newSig
        }
    }
    }

    return (sig_data)
}

#' Reads precomputed signature values form a tab-delimited text file
#'
#' First row of file contains sample labels that the signatures correspond with
#'
#' Each subsequent row contains a signature name in the first column,
#' followed by the signature type (either 'numerical' or 'factor')
#' followed by the signature values, one for each sample label in the file
#'
#' @importFrom utils read.table
#' @param filename Filename to read the precomputed signatures from
#' @param sampleLabels List of labels for which we want the signature scores
#' @param sep seperator to use when reading in file. Default is tab.
#' @return List of Signature data types
readPrecomputed <- function(filename, sampleLabels, sep="\t") {


    message("Loading data from ", filename)

    f <- as.matrix(read.table(filename, sep=sep))
    l1 <- f[1,]
    l2 <- f[2,]

    if ( (length(l1) != length(l2)) && (length(l1) != length(l2)-2)) {
    stop("Error in header line of precomputed signature file.\n
            First row should contain tab-separated list of samples")
    }

    ## Make sure that l1 is in the same order as sample_labels; match indices
    ## between signatrues in file & sampleLabels
    sampleLabels <- as.character(as.matrix(sampleLabels))

    common <- intersect(l1[-c(1,2)], sampleLabels)
    if (length(common) != length(sampleLabels)){
        stop("Provided precomputed signature dataframe must have same sample labels as the expression matrix")
    }

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
        stop("Error in sample ", sampleLabels[ind-2])
        })

    } else if (sigType == "factor") {
        sigIsFactor <- TRUE
        sigVals <- as.factor(as.matrix(sigVals))

    } else {
        stop("Column 2 of precomputed signature file should specify either 'numerical' or 'factor'")
    }

    sigScores <- c(sigScores, SignatureScores(sigVals, sigName, sampleLabels,
                                                sigIsFactor, TRUE, 0))
    }

    return(sigScores)

}

