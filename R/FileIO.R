#' Reads in a list of signature input files.
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

    if (!file.exists(filename)) {
        stop(paste0("Cannot find file: ", filename))
    }

    fsplit <- strsplit(basename(filename), "\\.")[[1]]
    file_ext <- fsplit[length(fsplit)]

    if (file_ext == "gmt") {

        conn <- file(filename, open="r")
        lines <-readLines(conn)
        close(conn)

        lines <- strsplit(lines, '\t')

        # Remove invalid lines
        lines <- lines[vapply(lines,
                              function(x) length(x) >= 3, TRUE)]

        for (line in lines) {

            genes <- c()
            values <- c()

            sig <- line[1]
            description <- line[2]
            gene_fields <- line[3:length(line)]

            sig_parts <- unlist(strsplit(sig, "_"))
            key <- tolower(sig_parts[[length(sig_parts)]])

            if (is.element(key, names(keys))) {
                sig_name <- paste(sig_parts[1:length(sig_parts)-1],
                                  collapse="_")
            } else {
                key <- "none"
                sig_name <- sig
            }


            for (k in gene_fields) {
                elem <- strsplit(k, ",")[[1]]

                if (length(elem) > 1) { # if the format is gene,value
                    gene <- elem[1]
                    value <- as.numeric(elem[2])
                } else {
                    gene <- elem[1]
                    value <- keys[key]
                }

                genes <- c(genes, gene)
                values <- c(values, value)
            }

            names(values) <- genes

            if (sig_name %in% names(sig_data)) {
                existingSig <- sig_data[[sig_name]]
                existingSig@sigDict[names(values)] <- values
                sig_data[[sig_name]] <- existingSig
            } else {
                newSig <- Signature(values, sig_name,
                                    filename, metaData=description)
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
        newSig <- Signature(sigList[[sig]], sig, filename, "")
        sig_data[[sig]] <- newSig
        }
    }
    }

    return (sig_data)
}
