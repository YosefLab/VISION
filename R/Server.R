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

#' Wrapper class for the P value matrix calculated during sigsVsProjections
#'
#' @param data P values for each signature, projection pair in the form of a matrix
#' @param proj_labels Projection names
#' @param sig_labels Names of signatures
#' @return ServerPMatrix object
ServerPMatrix <- function(data, proj_labels, sig_labels) {
            .Object <- new("ServerPMatrix", data=data, proj_labels=proj_labels,
                           sig_labels=sig_labels)
            return(.Object)
            }

#' Converts Signature object to JSON
#' @importFrom jsonlite toJSON
#' @param sig Signature object
#' @return JSON formatted Signature object.
signatureToJSON <- function(sig) {

    # Pass in a Signature object from a FastProject Object to be converted into a JSON object
    sig@sigDict <- as.list(sig@sigDict)

    json <- toJSON(sig, force=TRUE, pretty=TRUE, auto_unbox=TRUE)
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
sigScoresToJSON <- function(ss) {

    s <- as.list(ss[[1]])
    names(s) <- rownames(ss)
    json <- toJSON(s, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

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

#' Converts a sigProjMatrix from a FastProject Object to a JSON object
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

    json <- toJSON(sSPM, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

    return(json)
}

#' convert perason correlation coeffcients between PCs and sgnatures into a JSON object
#' @param pc the pearson correlation coefficients matrix
#' @param sigs the signatures of interest
#' @return Subsetted pearson correlations converted to JSON
pearsonCorrToJSON <- function(pc, sigs) {

    pc <- pc[sigs, ,drop = FALSE]
    cn <- paste("PC", 1:ncol(pc))

    sPC <- ServerSigProjMatrix(unname(pc), unname(pc), cn, sigs)

    json <- toJSON(sPC, force = TRUE, pretty = TRUE, auto_unbox = TRUE)

    return(json)

}

#' Convert a Cluster object to JSON
#' @importFrom jsonlite toJSON
#' @param cluster Cluster object
#' @return Cluster object converted to JSON
clusterToJSON <- function(cluster) {

    out <- list()
    out[['method']] <- cluster@method
    out[['param']] <- cluster@param
    out[['centers']] <- cluster@centers
    out[['data']] <- as.list(cluster@data[1,])
    json <- toJSON(out, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

    return(json)
}

#' Run the analysis again wth user-defined subsets or confguration
#' @param nfp the new FastProject object to analyze
#' @return None
newAnalysis <- function(nfp) {
    saveAndViewResults(analyze(nfp))
}

#' Lanch the server
#' @importFrom jsonlite fromJSON
#' @importFrom utils browseURL URLdecode stack
#' @importFrom SummarizedExperiment values
#' @param object FastProject object or path to a file containing such an
#' object (saved using saveAndViewResults, or directly using saveRDS)
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @return None
launchServer <- function(object, port=NULL, host=NULL, browser=TRUE) {

    if (is.null(port)) {
        port <- sample(8000:9999, 1)
    }
    if (is.null(host)) {
        host <- "127.0.0.1"
    }

    # Load the static file whitelist
    whitelist_file <- system.file("html_output/whitelist.txt",
                                  package = "FastProjectR")
    static_whitelist <- scan(whitelist_file, what = "",
                             quiet = TRUE)

    url <- paste0("http://", host, ":", port, "/Results.html")
    message(paste("Navigate to", url, "in a browser to interact with the app."))

    if(browser){
        browseURL(url)
    }

    # Launch the server
    jug() %>%
      get("/Signature/Scores/(?<sig_name1>.*)", function(req, res, err) {
        sigMatrix <- object@sigScores
        name <- URLdecode(req$params$sig_name1)
        out <- "Signature does not exist!"
        if (name %in% colnames(sigMatrix)) {
          ss <- sigMatrix[, name, drop = FALSE]
          ss <- as.data.frame(ss)
          out <- FastProjectR:::sigScoresToJSON(ss)
        }
        return(out)
      }) %>%
      get("/Signature/Meta/(?<sig_name3>.*)", function(req, res, err) {
        metaData <- object@metaData
        name <- URLdecode(req$params$sig_name3)
        out <- "Signature does not exist!"
        if (name %in% colnames(metaData)) {
          out <- FastProjectR:::sigScoresToJSON(metaData[name])
        }
        return(out)
      }) %>%
      get("/Signature/ListMeta", function(req, res, err){
        signatures <- object@sigData
        keys <- lapply(signatures, function(x) x@name)
        vals <- lapply(signatures, function(x) x@isMeta)
        names(vals) <- keys
        out <- toJSON(vals, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/Signature/Info/(?<sig_name2>.*)", function(req, res, err){
        signatures <- object@sigData
        name <- URLdecode(req$params$sig_name2)
        out <- "Signature does not exist!"
        if (name %in% names(signatures)) {
          sig <- signatures[[name]]
          out <- FastProjectR:::signatureToJSON(sig)
        }
        return(out)
      }) %>%
      get("/Signature/Expression/(?<sig_name4>.*)", function(req, res, err) {
        all_names = vapply(object@sigData, function(x){ return(x@name) }, "")
        name <- URLdecode(req$params$sig_name4)
        index = match(name, all_names)
        if(is.na(index)){
            out <- "Signature does not exist!"
        }
        else{
            sig = object@sigData[[index]]
            if(sig@isMeta) {
                stop("Can't get expression for meta data signature")
            }
            genes = names(sig@sigDict)
            expMat = object@exprData
            return(FastProjectR:::expressionToJSON(expMat, genes, zscore=TRUE))
        }
        return(out)
      }) %>%
      get("/FilterGroup/SigClusters/Normal", function(req, res, err) {

        cls <- object@ProjectionData@sigClusters
        cls <- cls$Computed

        out <- toJSON(cls, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/SigClusters/Meta", function(req, res, err) {

        cls <- object@ProjectionData@sigClusters
        cls <- cls$Meta

        out <- toJSON(cls, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/Projections/(?<proj_name1>.*)/coordinates", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name1)
        out <- FastProjectR:::coordinatesToJSON(object@Projections[[proj]])
        return(out)
      }) %>%
      get("/Projections/list", function(req, res, err) {
        proj_names <- names(object@Projections)
        out <- toJSON(proj_names, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/Projections/SigProjMatrix/Normal", function(req, res, err) {

        sigs <- colnames(object@sigScores)

        out <- FastProjectR:::sigProjMatrixToJSON(
                                  object@ProjectionData@sigProjMatrix,
                                  object@ProjectionData@emp_pMatrix,
                                  sigs)
        return(out)
      }) %>%
      get("/Projections/SigProjMatrix/Meta", function(req, res, err) {

        sigs <- colnames(object@metaData)

        out <- FastProjectR:::sigProjMatrixToJSON(
                                  object@ProjectionData@sigProjMatrix,
                                  object@ProjectionData@emp_pMatrix,
                                  sigs)
        return(out)
      }) %>%
      get("/Clusters/SigProjMatrix/Normal", function(req, res, err) {

        sigs <- colnames(object@sigScores)

        out <- FastProjectR:::sigProjMatrixToJSON(
                                  object@ClusterProjectionData@sigProjMatrix,
                                  object@ClusterProjectionData@emp_pMatrix,
                                  sigs)
        return(out)
      }) %>%
      get("/Clusters/SigProjMatrix/Meta", function(req, res, err) {

        sigs <- colnames(object@metaData)

        out <- FastProjectR:::sigProjMatrixToJSON(
                                  object@ClusterProjectionData@sigProjMatrix,
                                  object@ClusterProjectionData@emp_pMatrix,
                                  sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<proj_name2>.*)/clusters/(?<cluster_procedure>.*)/(?<param>.*)", function(req, res, err) {
        # projData <- object@projData

        proj <- URLdecode(req$params$proj_name2)
        method <- URLdecode(req$params$cluster_procedure)
        param <- as.numeric(URLdecode(req$params$param))

        clust <- FastProjectR:::cluster(object@ProjectionData@projections[[proj]], method, param)
        out <- FastProjectR:::clusterToJSON(clust)
        return(out)
      }) %>%
      get("/Tree/List", function(req, res, err) {

        ## all Trees have the same adjacency matrix, so we can use the first one
        W <- object@TreeProjectionData@projections[[1]]@adjMat

        return(toJSON(W))
      }) %>%
      get("/Tree/(?<proj_name3>.*)/Points", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name3)

        C <- object@TreeProjectionData@projections[[proj]]@vData

        return(toJSON(C))
      }) %>%
      get("/Tree/(?<proj_name4>.*)/Projection", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name4)

        out <- FastProjectR:::coordinatesToJSON(object@TreeProjectionData@projections[[proj]]@pData)

        return(out)
      }) %>%
      get("/Tree/SigProjMatrix/Normal", function(req, res, err) {

        sigs <- colnames(object@sigScores)
        out <- FastProjectR:::sigProjMatrixToJSON(
                                  object@TreeProjectionData@sigProjMatrix,
                                  object@TreeProjectionData@emp_pMatrix,
                                  sigs)
        return(out)
      }) %>%
      get("/Tree/SigProjMatrix/Meta", function(req, res, err) {

        sigs <- colnames(object@metaData)
        out <- FastProjectR:::sigProjMatrixToJSON(
                                  object@TreeProjectionData@sigProjMatrix,
                                  object@TreeProjectionData@emp_pMatrix,
                                  sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<pc_num1>.*)/Loadings/Positive", function(req, res, err) {
        pcnum <- as.numeric(URLdecode(req$params$pc_num1))

        c <- object@latentSpace[, pcnum, drop = FALSE] # cells x 1
        edata <- object@exprData # genes x cells
        l <- edata %*% c
        names(l) <- rownames(edata)

        l <- sqrt(sum(l ^ 2)) # normalize the vector

        posl <- l[l >= 0]

        js1 <- toJSON(with(stack(posl), tapply(values, ind, c, simplify = FALSE)))

        return(js1)

      }) %>%
      get("/FilterGroup/(?<pc_num6>.*)/Loadings/Negative", function(req, res, err) {
        pcnum <- as.numeric(URLdecode(req$params$pc_num6))

        c <- object@latentSpace[, pcnum, drop = FALSE] # cells x 1
        edata <- object@exprData # genes x cells
        l <- edata %*% c
        names(l) <- rownames(edata)

        l <- sqrt(sum(l ^ 2)) # normalize the vector

        negl <- l[l < 0]

        js2 <- toJSON(with(stack(negl),
                           tapply(values, ind, c, simplify = FALSE)))

        return(js2)

      }) %>%
      get("/FilterGroup/PCA/Coordinates", function(req, res, err) {

        pc <- object@latentSpace
        out <- FastProjectR:::coordinatesToJSON(pc)

        return(out)
      }) %>%
      get("/FilterGroup/PCVersus/(?<pc_num3>.*)/(?<pc_num4>.*)", function(req, res, err) {

        pc1 <- as.numeric(URLdecode(req$params$pc_num3))
        pc2 <- as.numeric(URLdecode(req$params$pc_num4))


        pcdata1 <- object@latentSpace[, pc1]
        pcdata2 <- object@latentSpace[, pc2]


        ret <- cbind(pcdata1, pcdata2)
        coord <- apply(unname(ret), 1, as.list)
        names(coord) <- rownames(ret)

        return(toJSON(coord, force = TRUE, auto_unbox = TRUE))
      }) %>%
      get("/FilterGroup/PearsonCorr/Normal", function(req, res, err) {

          sigs <- colnames(object@sigScores)

          pc <- object@PCAnnotatorData@pearsonCorr[, 1:10]

        return(FastProjectR:::pearsonCorrToJSON(pc, sigs))
      }) %>%
      get("/FilterGroup/PearsonCorr/Meta", function(req, res, err) {

          sigs <- colnames(object@metaData)
          numericMeta <- vapply(sigs,
                                function(x) is.numeric(object@metaData[[x]]),
                                FUN.VALUE = TRUE)
          sigs <- sigs[numericMeta]

          pc <- object@PCAnnotatorData@pearsonCorr[, 1:10]

        return(FastProjectR:::pearsonCorrToJSON(pc, sigs))
      }) %>%
      get("/FilterGroup/PearsonCorr/list", function(req, res, err) {

          pc <- object@PCAnnotatorData@pearsonCorr
          pcnames <- seq(ncol(pc))[1:10]
          result <- toJSON(
                           pcnames,
                           force=TRUE, pretty=TRUE, auto_unbox=TRUE
                           )
          return(result)
      }) %>%
      get("/Expression/Genes/List", function(req, res, err) {

        data <- object@exprData
        genes = rownames(data)

        result <- toJSON(
                         genes,
                         force=TRUE, pretty=TRUE, auto_unbox=TRUE
                         )

        return(result)

      }) %>%
      get("/Expression/Gene/(?<gene_name2>.*)", function(req, res, err) {

        data <- log2(object@exprData + 1)
        gene_name <- URLdecode(req$params$gene_name2)

        result <- toJSON(
                         as.list(data[gene_name,]),
                         force=TRUE, pretty=TRUE, auto_unbox=TRUE
                         )

        return(result)

      }) %>%
      get("/Clusters", function(req, res, err) {
        metaData <- object@metaData
        name <- object@cluster_variable
        out <- "No Clusters!"
        if (name %in% colnames(metaData)) {
          out <- FastProjectR:::sigScoresToJSON(metaData[name])
        }
        return(out)
      }) %>%
      get("/SessionInfo", function(req, res, err) {

        info <- list()

        if (.hasSlot(object, "name") && !is.null(object@name)) {
            info["name"] <- object@name
        } else {
            info["name"] <- ""
        }

        W <- object@TreeProjectionData
        hasTree <- !is.null(W)

        info["has_tree"] <- hasTree

        result <- toJSON(
                         info,
                         force = TRUE, pretty = TRUE, auto_unbox = TRUE
                         )

        return(result)

      }) %>%
      post("/Analysis/Run/", function(req, res, err) {
        subset <- fromJSON(req$body)
        subset <- subset[!is.na(subset)]

        if (length(object@pools) > 0) {
            clust <- object@pools[subset]
            subset <- unlist(clust)
        }

        nfp <- FastProjectR:::createNewFP(object, subset)
        FastProjectR:::newAnalysis(nfp)
        return()
      }) %>%
      get(path = NULL, function(req, res, err) {

          if (req$path == "/") {
              path <- "Results.html"
          } else {
              path <- substring(req$path, 2) # remove first / character
          }
          path <- paste0("html_output/", path)
          file_index <- match(path, static_whitelist)

          if (is.na(file_index)) {
              res$set_status(404)
              return(NULL)
          }

          file_path <- system.file(static_whitelist[file_index],
                         package = "FastProjectR")

          mime_type <- mime::guess_type(file_path)
          res$content_type(mime_type)

          data <- readBin(file_path, "raw", n = file.info(file_path)$size)

          if (grepl("image|octet|pdf", mime_type)) {
              return(data)
          } else {
              return(rawToChar(data))
          }
      }) %>%
      simple_error_handler_json() %>%
      serve_it(host=host, port=port)
}
