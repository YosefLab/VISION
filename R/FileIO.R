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
    file_ext <- fsplit[length(fsplit)]

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
            newSig <- Signature(values, sig_name, filename, "")
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
        newSig <- Signature(sigList[[sig]], sig, filename, "", isMeta=FALSE,
                            isFactor=FALSE)
        sig_data[[sig]] <- newSig
        }
    }
    }

    return (sig_data)
}
