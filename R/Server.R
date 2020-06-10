#' Wrapper class for gene expression object for JSON.
#'
#' @param data numeric matrix
#' @param sample_labels Labels of samples in expression data
#' @param gene_labels Lables of genes in expression data
#' @return ServerExpression object
ServerExpression <- function(data, sample_labels, gene_labels) {
            .Object <- new("ServerExpression", data=data,
                           sample_labels=sample_labels,
                           gene_labels=gene_labels)
            return(.Object)
            }

#' Wrapper class for Signature Projection Matrix
#'
#' @param zscores numeric matrix
#' @param pvals numeric matrix
#' @param proj_labels Projection names
#' @param sig_labels Names of signatures
#' @return ServerSigProjMatrix object
ServerSigProjMatrix <- function(zscores, pvals, proj_labels, sig_labels) {
            .Object <- new("ServerSigProjMatrix",
                           zscores = zscores, pvals = pvals,
                           proj_labels = proj_labels, sig_labels = sig_labels)
            return(.Object)
            }

#' Converts Signature object to JSON
#' @importFrom jsonlite toJSON
#' @param sig Signature object
#' @param geneImportance named numeric vector with gene correlations
#' @return JSON formatted Signature object.
signatureToJSON <- function(sig, geneImportance) {

    # Pass in a Signature object from an R Object to be converted into a JSON object
    out <- list()

    out$sigDict <- as.list(sig@sigDict)
    out$name <- sig@name
    out$source <- sig@source
    out$metaData <- sig@metaData
    out$geneImportance <- as.list(geneImportance)

    json <- toJSON(
        out, force = TRUE, pretty = TRUE,
        auto_unbox = TRUE, digits = I(4))
    return(json)

}

#' Converts row of sigantures score matrix to JSON
#' @importFrom jsonlite toJSON
#' @param names character vector of labels for signatures
#' @param values numeric vector for signature values
#' @return Signature scores list to JSON, with names of each entry that of the list names
sigScoresToJSON <- function(names, values) {

    out <- list(cells = names, values = values)

    json <- toJSON(
        out, force = TRUE, pretty = TRUE,
        auto_unbox = TRUE, digits = I(4)
    )

    return(json)
}

#' Converts a projection into a JSON object mapping each sample to a projection coordinate.
#' @importFrom jsonlite toJSON
#' @param p Projection coordinate data (NUM_SAMPLES x NUM_COMPONENTS)
#' @return JSON object mapping each sample to a projection coordinate.
coordinatesToJSON <- function(p) {

    coord <- as.data.frame(p)

    # This switching is needed because toJSON will drop row labels
    # if they are integers for some reason
    coord["sample_labels"] <- rownames(coord)
    rownames(coord) <- NULL

    json <- toJSON(
        coord, force = TRUE, pretty = TRUE,
        auto_unbox = TRUE, dataframe = "values",
        digits = I(4)
    )

    return(json)
}

#' Converts a sigProjMatrix from an R Object to a JSON object
#' @importFrom jsonlite toJSON
#' @param sigzscores Matrix of signature z-scores
#' @param sigpvals Matrix of signature p-values
#' @param sigs Signatures to subset zscores/pvalues
#' @return Subsetted sigProjMatrix converted to JSON
sigProjMatrixToJSON <- function(sigzscores, sigpvals, sigs) {

    sigs <- intersect(sigs, rownames(sigzscores))
    sigpvals <- sigpvals[sigs, , drop = FALSE]

    sigzscores <- sigzscores[rownames(sigpvals), , drop = FALSE]
    sigzscores <- sigzscores[, colnames(sigpvals), drop = FALSE]

    sSPM <- ServerSigProjMatrix(unname(sigzscores), unname(sigpvals), colnames(sigpvals), rownames(sigpvals))

    json <- toJSON(sSPM, force = TRUE, pretty = TRUE, digits = I(4))

    return(json)
}

#' convert perason correlation coeffcients between PCs and sgnatures into a JSON object
#' @param pc the pearson correlation coefficients matrix
#' @param sigs the signatures of interest
#' @return Subsetted pearson correlations converted to JSON
pearsonCorrToJSON <- function(pc, sigs) {

    pc <- pc[sigs, , drop = FALSE]
    cn <- colnames(pc)

    if (is.null(cn)){
        cn <- paste0("Comp ", seq_len(ncol(pc)))
    }

    sPC <- ServerSigProjMatrix(unname(pc), unname(pc), cn, sigs)

    json <- toJSON(sPC, force = TRUE, pretty = TRUE, digits = I(4))

    return(json)

}

compressJSONResponse <- function(req, res){

    res$headers[["Content-type"]] <- "application/json"

    if (!is.null(req$HTTP_ACCEPT_ENCODING) &&
        grepl("gzip", req$HTTP_ACCEPT_ENCODING, ignore.case = TRUE)
    ){
        res$setHeader("Content-Encoding", "gzip")
        GZIP_HEADER <- as.raw(c(31, 139, 8, 0, 0, 0, 0, 0, 4, 3))
        compressed <- memCompress(charToRaw(res$body))
        compressed <- compressed[-c(1, 2)]
        compressed <- compressed[- (
                (length(compressed) - 3):length(compressed)
        )]
        compressed <- c(GZIP_HEADER, compressed)
        res$body <- compressed
    }
}

#' Lanch the server
#' @importFrom jsonlite fromJSON
#' @importFrom utils browseURL URLdecode stack
#' @importFrom plumber plumber forward
#' @importFrom Matrix rowMeans
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats setNames
#' @importFrom mime guess_type
#' @param object Vision object
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @param browser Whether or not to launch the web browser
#' @return object
launchServer <- function(object, port=NULL, host=NULL, browser=TRUE) {

    if (is.null(port)) {
        port <- sample(8000:9999, 1)
    }
    if (is.null(host)) {
        host <- "127.0.0.1"
    }

    # Make sure all projections have column names
    projections <- object@Projections
    n <- names(projections)
    projections <- lapply(setNames(n, n), function(pname){
        proj <- projections[[pname]]
        if (is.null(colnames(proj))){
            colnames(proj) <- paste0(pname, "-", seq_len(ncol(proj)))
        }
        return(proj)
    })

    if (object@version < 1.11 && "Latent Space" %in% names(projections)){
        projections[["Latent Space"]] <- NULL
    }
    object@Projections <- projections

    # Make sure latent space columns have names
    if (is.null(colnames(object@LatentSpace))) {
        colnames(object@LatentSpace) <- paste0("Comp ", seq_len(ncol(object@LatentSpace)))
    }

    # Ensure the 'Viewer' list is fully initialized
    if (!("selections" %in% names(object@Viewer))){
        object@Viewer[["selections"]] <- list()
    }
    if (!("de_cache" %in% names(object@Viewer))){
        object@Viewer[["de_cache"]] <- list()
    }


    # Load the static file whitelist
    whitelist_file <- system.file("html_output/whitelist.txt",
                                  package = "VISION")
    static_whitelist <- scan(whitelist_file, what = "",
                             quiet = TRUE)

    url <- paste0("http://", host, ":", port, "/Results.html")
    message(paste("Navigate to", url, "in a browser to interact with the app."))

    if(browser){
        browseURL(url)
    }

    # Define the API

    pr <- plumber$new()

    pr$filter("defaultHeaders", function(req, res){
        res$setHeader("Cache-Control", "max-age=7200")
        res$setHeader("Content-type", "application/json")
        forward()
    })

    pr$handle("GET", "/Signature/Scores/<sig_name>",
        function(req, res, sig_name) {

        sigMatrix <- object@SigScores
        name <- URLdecode(sig_name)
        if (name %in% colnames(sigMatrix)) {
            values <- sigMatrix[, name]
            names <- rownames(sigMatrix)
            res$body <- sigScoresToJSON(names, values)
            compressJSONResponse(req, res)
            return(res)
        } else {
            return("Signature does not exist!")
        }
    })

    pr$handle("GET", "/Signature/Meta/<sig_name>",
        function(req, res, sig_name) {

      metaData <- object@metaData
      name <- URLdecode(sig_name)
      if (name %in% colnames(metaData)) {
          names <- rownames(metaData)
          values <- metaData[[name]]
          res$body <- sigScoresToJSON(names, values)
          compressJSONResponse(req, res)
          return(res)
      } else {
          return("Signature does not exist!")
      }
    })

    pr$handle("GET", "/Signature/Info/<sig_name>",
        function(req, res, sig_name){

      signatures <- object@sigData
      name <- URLdecode(sig_name)
      out <- "Signature does not exist!"
      if (name %in% names(signatures)) {
        sig <- signatures[[name]]

          if (.hasSlot(object, "SigGeneImportance")) {
              geneImportance <- object@SigGeneImportance[[name]]
          } else {
              geneImportance <- sig@sigDict * 0
          }

        out <- signatureToJSON(sig, geneImportance)
      }
      res$body <- out
      return(res)
    })

    pr$handle("GET", "/Signature/Expression/<sig_name>/<cluster_var>",
        function(req, res, sig_name, cluster_var) {

        all_names <- vapply(object@sigData, function(x) x@name, "")
        name <- URLdecode(sig_name)
        index <- match(name, all_names)
        if (is.na(index)){
            return("Signature does not exist!")
        }
        else{
            sig <- object@sigData[[index]]
            genes <- names(sig@sigDict)
            expMat <- object@exprData[genes, ]
            # transpose to aggregate
            # note only agg over genes in sig
            # on click gene plots gene
            expMat <- t(expMat)
            clusters <- object@metaData[rownames(expMat), cluster_var]
            clusters <- as.factor(clusters)
            onehot <- lapply(setNames(levels(clusters), levels(clusters)),
                function(level) ifelse(clusters == level, 1, 0)
            )
            onehot <- t(do.call(cbind, onehot))
            onehot <- onehot / rowSums(onehot)
            x <- onehot %*% expMat

            # get rid of aggregate artifact and redo transpose
            y <- as.matrix(t(x))
            rownames(y) <- colnames(expMat)

            # z-score before cluster
            rm <- matrixStats::rowMeans2(y)
            rsd <- matrixStats::rowSds(y)
            rsd[rsd == 0] <- 1.0
            y <- (y - rm) / rsd

            # Cluster the rows
            d <- dist(y, method = "euclidean")
            fit <- hclust(d, method = "average")
            y <- y[fit$order, , drop = FALSE]

            # Cluster the columns
            d <- dist(t(y), method = "euclidean")
            fit <- hclust(d, method = "average")
            y <- y[, fit$order, drop = FALSE]

            sExpr <- ServerExpression(y, colnames(y), rownames(y))
            json_out <- toJSON(
                sExpr, force = TRUE, pretty = TRUE, auto_unbox = TRUE,
                digits = I(4)
            )

            res$body <- json_out
            compressJSONResponse(req, res)

            return(res)
        }
    })

    pr$handle("GET", "/Proteins/<protein_name>/Values",
        function(req, res, protein_name) {

        proteinMatrix <- object@proteinData
        name <- URLdecode(protein_name)
        if (name %in% colnames(proteinMatrix)) {
            values <- proteinMatrix[, name]
            names <- rownames(proteinMatrix)
            res$body <- sigScoresToJSON(names, values)
            compressJSONResponse(req, res)
            return(res)
        } else {
            return("Proteins does not exist!")
        }
    })

    pr$handle("GET", "/FilterGroup/SigClusters/Normal", function(req, res) {

        cls <- object@LocalAutocorrelation$Clusters
        cls <- cls$Computed

        res$body <- toJSON(cls, auto_unbox = TRUE)
        return(res)

    })

    pr$handle("GET", "/FilterGroup/SigClusters/Proteins", function(req, res) {

        # Everything is in cluster 1
        cls <- as.list(numeric(nrow(object@LocalAutocorrelation$Proteins)) + 1)
        names(cls) <- rownames(object@LocalAutocorrelation$Proteins)

        res$body <- toJSON(cls, auto_unbox = TRUE)
        return(res)
    })

    pr$handle("GET", "/FilterGroup/SigClusters/Meta", function(req, res) {

        cls <- object@LocalAutocorrelation$Clusters
        cls <- cls$Meta

        res$body <- toJSON(cls, auto_unbox = TRUE)
        return(res)
    })

    pr$handle("GET", "/Projections/<proj_name>/coordinates/<proj_col>",
        function(req, res, proj_name, proj_col) {

        proj <- URLdecode(proj_name)
        col <- URLdecode(proj_col)

        if (proj %in% names(object@Projections)) {
            coords <- object@Projections[[proj]][, col, drop = FALSE]
        } else if (proj == "Meta Data") {
            coords <- object@metaData[, col, drop = FALSE]
        } else {
            coords <- object@LatentSpace[, col, drop = FALSE]
        }

        res$body <- coordinatesToJSON(coords)
        compressJSONResponse(req, res)

        return(res)

    })

    pr$handle("GET", "/Projections/list",
        function(req, res) {

        proj_names <- lapply(object@Projections, colnames)

        # Add in Latent Space
        latentSpaceName <- getParam(object, "latentSpace", "name")
        proj_names[[latentSpaceName]] <- colnames(object@LatentSpace)

        # Add in Meta-Data if it's numeric
        numericMeta <- vapply(seq_len(ncol(object@metaData)),
                              function(i) is.numeric(object@metaData[[i]]),
                              FUN.VALUE = TRUE)
        numericMeta <- colnames(object@metaData)[numericMeta]
        if (length(numericMeta) > 1){
            proj_names[["Meta Data"]] <- numericMeta
        }

        proj_names <- proj_names[order(names(proj_names), decreasing = TRUE)] # hack to make tsne on top
        res$body <- toJSON(proj_names)
        return(res)

    })

    pr$handle("GET", "/Tree/Projections/list",
        function(req, res) {

        if (is.null(object@TrajectoryProjections)){
            proj_names <- character()
        } else {
            proj_names <- names(object@TrajectoryProjections)
        }

        res$body <- toJSON(proj_names)
        return(res)
    })

    pr$handle("GET", "/Tree/Projections/<proj_name>/coordinates",
        function(req, res, proj_name) {

        proj <- URLdecode(proj_name)

        coords <- object@TrajectoryProjections[[proj]]@pData
        C <- object@TrajectoryProjections[[proj]]@vData
        W <- object@TrajectoryProjections[[proj]]@adjMat

        coords <- as.data.frame(coords)
        coords["sample_labels"] <- rownames(coords)
        rownames(coords) <- NULL

        out <- list(coords, C, W)

        res$body <- toJSON(
            out, force = TRUE, pretty = TRUE,
            auto_unbox = TRUE, dataframe = "values", digits = I(4)
        )

        compressJSONResponse(req, res)

        return(res)

    })

    pr$handle("GET", "/Tree/SigProjMatrix/Normal", function(req, res) {

        metaData <- object@metaData
        meta_n <- vapply(names(metaData),
            function(metaName) is.numeric(metaData[, metaName]),
            FUN.VALUE = TRUE)

        meta_n <- colnames(metaData)[meta_n]

        data <- object@TrajectoryAutocorrelation$Signatures
        stat_s <- as.matrix(data[, "C", drop = F])
        pval_s <- as.matrix(data[, "FDR", drop = F])

        data <- object@TrajectoryAutocorrelation$Meta[meta_n, , drop = F]
        stat_p <- as.matrix(data[, "C", drop = F])
        pval_p <- as.matrix(data[, "FDR", drop = F])

        stat <- rbind(stat_s, stat_p)
        pval <- rbind(pval_s, pval_p)

        colnames(stat) <- "Score"
        colnames(pval) <- "Score"

        res$body <- sigProjMatrixToJSON(stat, pval, sigs = rownames(stat))

        return(res)
    })

    pr$handle("GET", "/Tree/ProteinMatrix", function(req, res) {

        data <- object@TrajectoryAutocorrelation$Proteins
        stat <- as.matrix(data[, "C", drop = F])
        pval <- as.matrix(data[, "FDR", drop = F])

        colnames(stat) <- "Score"
        colnames(pval) <- "Score"

        res$body <- sigProjMatrixToJSON(stat, pval, sigs = rownames(stat))

        return(res)
    })

    pr$handle("GET", "/Tree/SigProjMatrix/Meta", function(req, res) {

        data <- object@TrajectoryAutocorrelation$Meta
        stat <- as.matrix(data[, "C", drop = F])
        pval <- as.matrix(data[, "FDR", drop = F])

        colnames(stat) <- "Score"
        colnames(pval) <- "Score"

        res$body <- sigProjMatrixToJSON(stat, pval, sigs = rownames(stat))

        return(res)
    })

    pr$handle("GET", "/PearsonCorr/Normal", function(req, res) {

        pc <- object@LCAnnotatorData@pearsonCorr
        sigs <- rownames(pc)

        res$body <- pearsonCorrToJSON(pc, sigs)

        return(res)
    })

    pr$handle("GET", "/PearsonCorr/Proteins", function(req, res) {

        pc <- object@LCAnnotatorData@pearsonCorrProteins
        sigs <- rownames(pc)

        res$body <- pearsonCorrToJSON(pc, sigs)

        return(res)
    })

    pr$handle("GET", "/PearsonCorr/Meta", function(req, res) {

        sigs <- colnames(object@metaData)
        numericMeta <- vapply(sigs,
                              function(x) is.numeric(object@metaData[[x]]),
                              FUN.VALUE = TRUE)
        sigs <- sigs[numericMeta]

        pc <- object@LCAnnotatorData@pearsonCorr

        res$body <- pearsonCorrToJSON(pc, sigs)
        return(res)
    })

    pr$handle("GET", "/Expression/Genes/List", function(req, res) {

        if (hasUnnormalizedData(object)) {
            data <- object@unnormalizedData
        } else {
            data <- object@exprData
        }
        genes <- rownames(data)

        res$body <- toJSON(genes, force = TRUE, pretty = TRUE)

        compressJSONResponse(req, res)

        return(res)
    })

    pr$handle("GET", "/Expression/Gene/<gene_name>",
        function(req, res, gene_name) {

        gene_name <- URLdecode(gene_name)

        if (hasUnnormalizedData(object)) {
            data <- object@unnormalizedData
        } else {
            data <- object@exprData
        }

        cells <- colnames(data)
        data <- log2(as.numeric(data[gene_name, ]) + 1)

        out <- list(cells = cells, values = data)

        res$body <- toJSON(
            out, force = TRUE, pretty = TRUE, auto_unbox = TRUE,
            digits = I(4)
        )

        compressJSONResponse(req, res)

        return(res)
    })

    pr$handle("GET", "/Clusters/list",
        function(req, res) {

        cluster_vars <- names(object@ClusterComparisons[[1]])

        res$body <- toJSON(cluster_vars, force = TRUE, pretty = TRUE)
        return(res)
    })
    
    pr$handle("POST", "/DE", function(req, res) {
        # Params
        body <- fromJSON(req$postBody)

        exprData <- object@exprData

        type_n <- body$type_n
        type_d <- body$type_d

        subtype_n <- body$subtype_n
        subtype_d <- body$subtype_d

        group_num <- body$group_num
        group_denom <- body$group_denom

        cells_min <- as.numeric(body$min_cells) # default 1
        subsample_groups <- body$subsample_groups
        subsample_cells_n <- as.numeric(body$subsample_N)

        # I want to skip caching a manual selection because the hash would be complicated
        skip_cache <- FALSE;
        if (type_n == "current") {
          skip_cache <- TRUE;
        } else {
          # hash by the params as a string
          hashed_body <- paste(type_n, type_d, subtype_n, subtype_d, group_num,
              group_denom, subsample_cells_n, cells_min, subsample_groups, sep = " ")
        }
        if (skip_cache || is.null(object@Viewer$de_cache[[hashed_body]])) {
            # The case where we know the de calculation is not cached, so we have to actually
            # calculate it.

            selections <- getSelections(object)

            if (type_n == "current") {
                cells_num <- unlist(strsplit(group_num, ","))
            } else if (type_n == "saved_selection") {
                cells_num <- selections[[group_num]]
            } else if (type_n == "meta") {
                cells_num <- rownames(object@metaData)[
                    which(object@metaData[[subtype_n]] == group_num)
                    ]
            } else {
                print("ERROR! Num type unrecognized: " + type_n)
            }

            if (type_d == "remainder") {
                cells_denom <- setdiff(colnames(exprData), cells_num)
            } else if (type_d == "saved_selection") {
                cells_denom <- selections[[group_denom]]
            } else if (type_d == "meta") {
                cells_denom <- rownames(object@metaData)[
                    which(object@metaData[[subtype_d]] == group_denom)
                    ]
            } else {
                print("ERROR! Denom type unrecognized: " + type_d)
            }

            # subsample the numerator and denominator groups
            if (subsample_groups) {
                if (length(cells_num) > subsample_cells_n) {
                    # if possible, the numerator becomes a random sample of size subsample_cells_n
                    cells_num <- sample(cells_num, size = subsample_cells_n)
                }
                if (length(cells_denom) > subsample_cells_n) {
                    # if possible, the denominator becomes a random sample of size subsample_cells_n
                    cells_denom <- sample(cells_denom, size = subsample_cells_n)
                }
            }

            # subset to only genes expressed in at least min cells
            if (cells_min < length(colnames(exprData)) && cells_min > 0) {
                exprDataSubset <- exprData[rowSums(exprData > 0) >= cells_min, , drop = FALSE]
            } else {
                exprDataSubset <- exprData
            }


            cluster_num <- match(cells_num, colnames(exprDataSubset))
            cluster_denom <- match(cells_denom, colnames(exprDataSubset))

            # Needed for sparse subsetting
            cluster_num <- as.numeric(cluster_num)
            cluster_denom <- as.numeric(cluster_denom)

            out <- matrix_wilcox_cpp(exprDataSubset, cluster_num, cluster_denom)

            out$stat <- pmax(out$AUC, 1 - out$AUC)
            out$Feature <- rownames(out)
            rownames(out) <- NULL

            numMean <- rowMeans(exprDataSubset[, cluster_num])
            denomMean <- rowMeans(exprDataSubset[, cluster_denom])
            bias <- 1 / sqrt(length(cluster_num) * length(cluster_denom))

            out$logFC <- log2( (numMean + bias) / (denomMean + bias) )

            out <- out[, c("Feature", "logFC", "stat", "pval"), drop = FALSE]
            out$Type <- "Gene"

            if (hasProteinData(object)) {
                proteinData <- t(object@proteinData)[, colnames(exprData), drop = FALSE]
                out_p <- matrix_wilcox_cpp(proteinData, cluster_num, cluster_denom)

                out_p$pval <- p.adjust(out_p$pval, method = "fdr")
                out_p$stat <- pmax(out_p$AUC, 1 - out_p$AUC)
                out_p$Feature <- rownames(out_p)
                rownames(out_p) <- NULL

                numMean <- rowMeans(proteinData[, cluster_num])
                denomMean <- rowMeans(proteinData[, cluster_denom])
                bias <- 1 / sqrt(length(cluster_num) * length(cluster_denom))

                out_p$logFC <- log2( (numMean + bias) / (denomMean + bias) )

                out_p <- out_p[, c("Feature", "logFC", "stat", "pval"), drop = FALSE]
                out_p$Type <- "Protein"

                out <- rbind(out, out_p)
            }

            out$pval <- p.adjust(out$pval, method = "fdr")

            # add the calculated results to the object, unless it was a user selection
            if (!skip_cache) {
                object@Viewer$de_cache[[hashed_body]] <<- out;
                # add to object body, removing the first stored item
                # TODO could implement a lru cache?
                # storing 20 max TODO decide on max size for cache
                if (length(object@Viewer$de_cache) > 20) {
                    # remove first cached item
                    object@Viewer$de_cache[1] <<- NULL;
                }
            }

        } else {
            # we have a cached result for these exact params, return it
            out <- object@Viewer$de_cache[[hashed_body]]
        }

        out <- as.list(out)

        result <- toJSON(
            out, force = TRUE, pretty = TRUE,
            auto_unbox = TRUE, use_signif = TRUE
        )
        res$body <- result
        compressJSONResponse(req, res)

        return(res)
    })

    pr$handle("GET", "/Clusters/<cluster_variable>/Cells",
        function(req, res, cluster_variable) {

        metaData <- object@metaData
        cluster_variable <- URLdecode(cluster_variable)
        if (cluster_variable %in% colnames(metaData)) {
            res$body <- sigScoresToJSON(
                names = rownames(metaData),
                values = metaData[[cluster_variable]]
                )
            compressJSONResponse(req, res)
            return(res)
        } else {
            return("No Clusters!")
        }
    })

    pr$handle("GET", "/Clusters/MetaLevels", function(req, res) {
        out <- list()
        metaData <- object@metaData
        for (cluster_variable in names(object@ClusterComparisons[[1]])) {
            out[[cluster_variable]] <- levels(metaData[[cluster_variable]])
        }

        res$body <- toJSON(
            out,
            force = TRUE, pretty = TRUE, auto_unbox = TRUE
        )
        compressJSONResponse(req, res)
        return(res)
    })

    pr$handle("GET", "/Clusters/<cluster_variable>/SigProjMatrix/Normal",
        function(req, res, cluster_variable) {

        # Gather signatures and numeric meta-data names
        sigs <- colnames(object@SigScores)

        metaData <- object@metaData
        meta_n <- vapply(names(metaData),
            function(metaName) is.numeric(metaData[, metaName]),
            FUN.VALUE = TRUE)

        meta_n <- colnames(metaData)[meta_n]

        pvals_s <- object@LocalAutocorrelation$Signatures[sigs, "FDR", drop = F]
        stat_s <- object@LocalAutocorrelation$Signatures[sigs, "C", drop = F]

        pvals_p <- object@LocalAutocorrelation$Meta[meta_n, "FDR", drop = F]
        stat_p <- object@LocalAutocorrelation$Meta[meta_n, "C", drop = F]

        pvals <- rbind(pvals_s, pvals_p)
        stat <- rbind(stat_s, stat_p)
        colnames(stat) <- "Score"
        colnames(pvals) <- "Score"
        stat <- as.matrix(stat)
        pvals <- as.matrix(pvals)

        cluster_variable <- URLdecode(cluster_variable)

        # Get Signature Cluster Comparisons
        var_res <- object@ClusterComparisons$Signatures[[cluster_variable]]

        cluster_pval_s <- lapply(var_res, function(var_level_res){
            var_level_res["FDR"]
        })
        cluster_pval_s <- as.matrix(do.call(cbind, cluster_pval_s))
        colnames(cluster_pval_s) <- names(var_res)

        cluster_stat_s <- lapply(var_res, function(var_level_res){
            var_level_res["stat"]
        })
        cluster_stat_s <- as.matrix(do.call(cbind, cluster_stat_s))
        colnames(cluster_stat_s) <- names(var_res)

        # Get Meta Cluster Comparisons
        var_res <- object@ClusterComparisons$Meta[[cluster_variable]]

        cluster_pval_p <- lapply(var_res, function(var_level_res){
            var_level_res["FDR"]
        })
        cluster_pval_p <- as.matrix(do.call(cbind, cluster_pval_p))
        colnames(cluster_pval_p) <- names(var_res)

        cluster_stat_p <- lapply(var_res, function(var_level_res){
            var_level_res["stat"]
        })
        cluster_stat_p <- as.matrix(do.call(cbind, cluster_stat_p))
        colnames(cluster_stat_p) <- names(var_res)

        # Combine together
        cluster_pval <- rbind(cluster_pval_s, cluster_pval_p)
        cluster_stat <- rbind(cluster_stat_s, cluster_stat_p)
        cluster_pval <- cluster_pval[rownames(pvals), , drop = F]
        cluster_stat <- cluster_stat[rownames(pvals), , drop = F]

        pvals <- cbind(pvals, cluster_pval)
        stat <- cbind(stat, cluster_stat)

        res$body <- sigProjMatrixToJSON(stat, pvals, sigs = rownames(stat))
        return(res)
    })

    pr$handle("GET", "/Clusters/<cluster_variable>/SigProjMatrix/Meta",
        function(req, res, cluster_variable) {

        cluster_variable <- URLdecode(cluster_variable)

        pvals <- object@LocalAutocorrelation$Meta[, "FDR", drop = F]
        stat <- object@LocalAutocorrelation$Meta[, "C", drop = F]
        colnames(stat) <- "Score"
        colnames(pvals) <- "Score"
        stat <- as.matrix(stat)
        pvals <- as.matrix(pvals)

        var_res <- object@ClusterComparisons$Meta[[cluster_variable]]

        cluster_pval <- lapply(var_res, function(var_level_res){
            var_level_res["FDR"]
        })
        cluster_pval <- as.matrix(do.call(cbind, cluster_pval))
        colnames(cluster_pval) <- names(var_res)

        cluster_stat <- lapply(var_res, function(var_level_res){
            var_level_res["stat"]
        })
        cluster_stat <- as.matrix(do.call(cbind, cluster_stat))
        colnames(cluster_stat) <- names(var_res)

        cluster_pval <- cluster_pval[rownames(pvals), , drop = F]
        cluster_stat <- cluster_stat[rownames(pvals), , drop = F]

        pvals <- cbind(pvals, cluster_pval)
        stat <- cbind(stat, cluster_stat)

        res$body <- sigProjMatrixToJSON(stat, pvals, sigs = rownames(stat))

        return(res)
    })

    pr$handle("GET", "/Clusters/<cluster_variable>/ProteinMatrix",
        function(req, res, cluster_variable) {

        proteins <- colnames(object@proteinData)

        cluster_variable <- URLdecode(cluster_variable)
        pvals <- object@LocalAutocorrelation$Proteins[, "FDR", drop = F]
        stat <- object@LocalAutocorrelation$Proteins[, "C", drop = F]
        colnames(stat) <- "Score"
        colnames(pvals) <- "Score"
        stat <- as.matrix(stat)
        pvals <- as.matrix(pvals)

        var_res <- object@ClusterComparisons$Proteins[[cluster_variable]]

        cluster_pval <- lapply(var_res, function(var_level_res){
            var_level_res["FDR"]
        })
        cluster_pval <- as.matrix(do.call(cbind, cluster_pval))
        colnames(cluster_pval) <- names(var_res)

        cluster_stat <- lapply(var_res, function(var_level_res){
            var_level_res["stat"]
        })
        cluster_stat <- as.matrix(do.call(cbind, cluster_stat))
        colnames(cluster_stat) <- names(var_res)

        cluster_pval <- cluster_pval[rownames(pvals), , drop = F]
        cluster_stat <- cluster_stat[rownames(pvals), , drop = F]

        pvals <- cbind(pvals, cluster_pval)
        stat <- cbind(stat, cluster_stat)

        res$body <- sigProjMatrixToJSON(stat, pvals, proteins)
        return(res)
    })

    pr$handle("GET", "/SessionInfo", function(req, res) {

        info <- list()

        info["name"] <- getParam(object, "name")

        W <- object@LatentTrajectory
        hasTree <- !is.null(W)

        info["has_tree"] <- hasTree

        info[["meta_sigs"]] <- colnames(object@metaData)

        info[["pooled"]] <- object@params$micropooling$pool

        info[["ncells"]] <- nrow(object@metaData)

        info[["has_sigs"]] <- length(object@sigData) > 0

        info[["has_proteins"]] <- hasProteinData(object)

        info[["has_lca"]] <- !is.null(object@LCAnnotatorData)
        
        if (!is.null(object@Tree)) {
          info[["dendrogram"]] <- write.tree(object@Tree)
        }
        
        info[["has_dendrogram"]] <- !is.null(object@Tree)

        res$body <- toJSON(info, force = TRUE,
            pretty = TRUE, auto_unbox = TRUE)

        return(res)

    })

    pr$handle("GET", "/Cell/<cell_id>/Meta",
        function(req, res, cell_id) {

        cell_id <- URLdecode(cell_id)
        cell_meta <- as.list(object@metaData[cell_id, ])
        res$body <- toJSON(cell_meta, auto_unbox = TRUE)

        return(res)
    })

    pr$handle("POST", "/Cells/Meta", function(req, res) {

        subset <- fromJSON(req$postBody)
        subset <- subset[!is.na(subset)]

        if (object@params$micropooling$pool){
            cells <- unname(unlist(object@Pools[subset]))
        } else {
            cells <- subset
        }

        metaSubset <- object@metaData[cells, ]

        numericMeta <- vapply(metaSubset, is.numeric, FUN.VALUE = TRUE)

        metaSummaryNumeric <-
            lapply(metaSubset[numericMeta], function(item){
               vals <- quantile(item, probs = c(0, .5, 1), na.rm = TRUE)
               names(vals) <- c("Min", "Median", "Max")
               return(as.list(vals))
            })

        metaSummaryFactor <-
            lapply(metaSubset[!numericMeta], function(item){
               vals <- sort(table(droplevels(item)),
                            decreasing = TRUE)
               vals <- vals / length(item) * 100
               return(as.list(vals))
            })

        metaSummary <- list(numeric = metaSummaryNumeric,
                            factor = metaSummaryFactor
                           )

        res$body <- toJSON(
            metaSummary, force = TRUE, pretty = TRUE,
            auto_unbox = TRUE, digits = I(4)
        )

        return(res)

    })

    pr$handle("GET", "/Cells/Selections", function(req, res) {

        # Disable caching for this request
        res$headers[["Cache-Control"]] <- NULL

        selections <- getSelections(object)

        selectionNames <- as.character(names(selections))
        res$body <- toJSON(selectionNames)
        return(res)
    })

    pr$handle("GET", "/Cells/Selections/<selection_id>",
        function(req, res, selection_id) {

        # Disable caching for this request
        res$headers[["Cache-Control"]] <- NULL

        selection_id <- URLdecode(selection_id)
        selections <- getSelections(object)
        selectionCells <- selections[[selection_id]]
        res$body <- toJSON(selectionCells, auto_unbox = TRUE)
        return(res)

    })

    pr$handle("POST", "/Cells/Selections/<selection_id>",
        function(req, res, selection_id) {

        selection_id <- URLdecode(selection_id)
        selections <- getSelections(object)
        cell_ids <- fromJSON(req$postBody)
        object@Viewer$selections[[selection_id]] <<- cell_ids # Super assignment!
        res$body <- "null"  # Empty body for successful POST
        return(res)
    })

    # Assume all other paths are files or 404

    pr$filter("filterFilter", function(req, res) {

        if (req$PATH_INFO == "/") {
            path <- "Results.html"
        } else {
            path <- substring(req$PATH_INFO, 2) # remove first / character
        }
        path <- paste0("html_output/", path)
        file_index <- match(path, static_whitelist)

        if (is.na(file_index)) {
            forward()
            return()
        }

        file_path <- system.file(static_whitelist[file_index],
                       package = "VISION")

        mime_type <- mime::guess_type(file_path)
        res$headers[["Content-type"]] <- mime_type

        data <- readBin(file_path, "raw", n = file.info(file_path)$size)

        if (grepl("image|octet|pdf", mime_type)) {
            res$body <- data
        } else {
            res$body <- rawToChar(data)
        }
        return(res)
    })

    tryCatch({
        pr$run(host = host, port = port, swagger = FALSE)
    },
    interrupt = function(i){
        message("Server Exited")
    })

    return(object)
}
