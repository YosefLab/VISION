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

    json <- toJSON(out, force=TRUE, pretty=TRUE, auto_unbox=TRUE)
    return(json)

}

#' Convertes expression matrix to JSON
#' @importFrom jsonlite toJSON
#' @param expr Expression Matrix
#' @param geneList optional list of genes to subset from expr
#' @return (Potentially) subsetted expression matrix
expressionToJSON <- function(expr, geneList=NULL, zscore=FALSE) {

    if (!is.null(geneList)) {
        geneList = intersect(geneList, rownames(expr))
        expr <- expr[geneList,, drop=FALSE]
    }

    expr <- as.matrix(expr)

    if (zscore) {
        rm <- matrixStats::rowMeans2(expr)
        rsd <- matrixStats::rowSds(expr)
        rsd[rsd == 0] <- 1.0
        expr_norm <- (expr - rm) / rsd
    }

    sExpr <- ServerExpression(expr_norm, colnames(expr), rownames(expr))

    ejson <- toJSON(sExpr, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

    return(ejson)
}

#' Converts row of sigantures score matrix to JSON
#' @importFrom jsonlite toJSON
#' @param ss single-column dataframe with scores for a single signature
#' @return Signature scores list to JSON, with names of each entry that of the list names
sigScoresToJSON <- function(names, values) {

    out <- list(cells = names, values = values)

    json <- toJSON(
                     out,
                     force = TRUE, pretty = TRUE, auto_unbox = TRUE
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

    json <- toJSON(coord, force=TRUE, pretty=TRUE,
                   auto_unbox=TRUE, dataframe="values")

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

    json <- toJSON(sSPM, force=TRUE, pretty=TRUE)

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

    json <- toJSON(sPC, force = TRUE, pretty = TRUE)

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
    if (is.null(colnames(object@latentSpace))) {
        colnames(object@latentSpace) <- paste0("Comp ", seq_len(ncol(object@latentSpace)))
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

        sigMatrix <- object@sigScores
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

    pr$handle("GET", "/Signature/Expression/<sig_name>",
        function(req, res, sig_name) {

        all_names <- vapply(object@sigData, function(x) x@name, "")
        name <- URLdecode(sig_name)
        index <- match(name, all_names)
        if (is.na(index)){
            return("Signature does not exist!")
        }
        else{
            sig <- object@sigData[[index]]
            genes <- names(sig@sigDict)
            expMat <- object@exprData
            res$body <- expressionToJSON(expMat, genes, zscore = TRUE)
            compressJSONResponse(req, res)
            return(res)
        }
    })

    pr$handle("GET", "/FilterGroup/SigClusters/Normal", function(req, res) {

        cls <- object@SigConsistencyScores@sigClusters
        cls <- cls$Computed

        res$body <- toJSON(cls, auto_unbox = TRUE)
        return(res)

    })

    pr$handle("GET", "/FilterGroup/SigClusters/Meta", function(req, res) {

        cls <- object@SigConsistencyScores@sigClusters
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
        } else {
            coords <- object@latentSpace[, col, drop = FALSE]
        }

        res$body <- coordinatesToJSON(coords)
        compressJSONResponse(req, res)

        return(res)

    })

    pr$handle("GET", "/Projections/list",
        function(req, res) {

        proj_names <- lapply(object@Projections, colnames)

        latentSpaceName <- getParam(object, "latentSpaceName")
        proj_names[[latentSpaceName]] <- colnames(object@latentSpace)

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

        res$body <- toJSON(out, force = TRUE, pretty = TRUE,
                           auto_unbox = TRUE, dataframe = "values")

        compressJSONResponse(req, res)

        return(res)

    })

    pr$handle("GET", "/Tree/SigProjMatrix/Normal", function(req, res) {

        sigs <- colnames(object@sigScores)

        metaData <- object@metaData
        meta_n <- vapply(names(metaData),
            function(metaName) is.numeric(metaData[, metaName]),
            FUN.VALUE = TRUE)

        meta_n <- colnames(metaData)[meta_n]
        sigs <- c(sigs, meta_n)

        res$body <- sigProjMatrixToJSON(
            object@TrajectoryConsistencyScores@Consistency,
            object@TrajectoryConsistencyScores@FDR,
            sigs)

        return(res)
    })

    pr$handle("GET", "/Tree/SigProjMatrix/Meta", function(req, res) {

        sigs <- colnames(object@metaData)
        res$body <- sigProjMatrixToJSON(
            object@TrajectoryConsistencyScores@Consistency,
            object@TrajectoryConsistencyScores@FDR,
            sigs)

        return(res)
    })

    pr$handle("GET", "/PearsonCorr/Normal", function(req, res) {

        pc <- object@PCAnnotatorData@pearsonCorr
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

        pc <- object@PCAnnotatorData@pearsonCorr

        res$body <- pearsonCorrToJSON(pc, sigs)
        return(res)
    })

    pr$handle("GET", "/PearsonCorr/list", function(req, res) {

        pc <- object@PCAnnotatorData@pearsonCorr
        pcnames <- seq(ncol(pc))[1:10]
        res$body <- toJSON(pcnames, force = TRUE, pretty = TRUE)
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

        data <- log2(as.numeric(data[gene_name, ]) + 1)

        out <- list(cells = names(data), values = data)

        res$body <- toJSON(
            out, force = TRUE, pretty = TRUE, auto_unbox = TRUE)

        compressJSONResponse(req, res)

        return(res)
    })

    pr$handle("GET", "/Clusters/list",
        function(req, res) {

        cluster_vars <- names(object@ClusterSigScores)
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

        if (type_n == "current") {
            cells_num <- unlist(strsplit(group_num, ","))
        } else if (type_n == "saved_selection") {
            cells_num <- object@selections[[group_num]]
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
            cells_denom <- object@selections[[group_denom]]
        } else if (type_d == "meta") {
            cells_denom <- rownames(object@metaData)[
                which(object@metaData[[subtype_d]] == group_denom)
                ]
        } else {
            print("ERROR! Denom type unrecognized: " + type_d)
        }


        cluster_num <- match(cells_num, colnames(exprData))
        cluster_denom <- match(cells_denom, colnames(exprData))

        out <- matrix_wilcox_cpp(exprData, cluster_num, cluster_denom)

        out$pval <- p.adjust(out$pval, method = "fdr")
        out$stat <- pmax(out$AUC, 1 - out$AUC)

        numMean <- rowMeans(exprData[, cluster_num])
        denomMean <- rowMeans(exprData[, cluster_denom])
        bias <- 1 / sqrt(length(cluster_num) * length(cluster_denom))

        out$logFC <- log2( (numMean + bias) / (denomMean + bias) )

        out <- out[, c("gene", "logFC", "stat", "pval"), drop = FALSE]
        out <- as.list(out)

        result <- toJSON(
          out,
          force = TRUE, pretty = TRUE, auto_unbox = TRUE, use_signif = TRUE
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
        for (cluster_variable in names(object@ClusterSigScores)) {
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

        sigs <- colnames(object@sigScores)

        metaData <- object@metaData
        meta_n <- vapply(names(metaData),
            function(metaName) is.numeric(metaData[, metaName]),
            FUN.VALUE = TRUE)

        meta_n <- colnames(metaData)[meta_n]
        sigs <- c(sigs, meta_n)

        cluster_variable <- URLdecode(cluster_variable)
        pvals <- object@SigConsistencyScores@FDR
        stat <- object@SigConsistencyScores@Consistency

        var_res <- object@ClusterSigScores[[cluster_variable]]

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

        res$body <- sigProjMatrixToJSON(stat, pvals, sigs)
        return(res)
    })

    pr$handle("GET", "/Clusters/<cluster_variable>/SigProjMatrix/Meta",
        function(req, res, cluster_variable) {

        sigs <- colnames(object@metaData)

        cluster_variable <- URLdecode(cluster_variable)
        pvals <- object@SigConsistencyScores@FDR
        stat <- object@SigConsistencyScores@Consistency

        var_res <- object@ClusterSigScores[[cluster_variable]]

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

        res$body <- sigProjMatrixToJSON(stat, pvals, sigs)

        return(res)
    })

    pr$handle("GET", "/SessionInfo", function(req, res) {

        info <- list()

        if (.hasSlot(object, "name") && !is.null(object@name)) {
            info["name"] <- object@name
        } else {
            info["name"] <- ""
        }

        W <- object@latentTrajectory
        hasTree <- !is.null(W)

        info["has_tree"] <- hasTree

        info[["meta_sigs"]] <- colnames(object@metaData)

        info[["pooled"]] <- object@pool

        info[["ncells"]] <- nrow(object@metaData)

        info[["has_sigs"]] <- length(object@sigData) > 0

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

        if (object@pool){
            cells <- unname(unlist(object@pools[subset]))
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
            metaSummary, force = TRUE, pretty = TRUE, auto_unbox = TRUE)

        return(res)

    })

    pr$handle("GET", "/Cells/Selections", function(req, res) {

        # Disable caching for this request
        res$headers[["Cache-Control"]] <- NULL

        selectionNames <- as.character(names(object@selections))
        res$body <- toJSON(selectionNames)
        return(res)
    })

    pr$handle("GET", "/Cells/Selections/<selection_id>",
        function(req, res, selection_id) {

        # Disable caching for this request
        res$headers[["Cache-Control"]] <- NULL

        selection_id <- URLdecode(selection_id)
        selectionCells <- object@selections[[selection_id]]
        res$body <- toJSON(selectionCells, auto_unbox = TRUE)
        return(res)

    })

    pr$handle("POST", "/Cells/Selections/<selection_id>",
        function(req, res, selection_id) {

        selection_id <- URLdecode(selection_id)
        cell_ids <- fromJSON(req$postBody)
        object@selections[[selection_id]] <<- cell_ids # Super assignment!
        res$body <- ""  # Empty body for successful POST
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
