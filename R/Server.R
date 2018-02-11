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
#' @param data Signatrue Projection matrix as obtained from sigsVsProjections
#' @param proj_labels Projection names
#' @param sig_labels Names of signatures
#' @return ServerSigProjMatrix object
ServerSigProjMatrix <- function(data, proj_labels, sig_labels) {
            .Object <- new("ServerSigProjMatrix", data=data,
                           proj_labels=proj_labels, sig_labels=sig_labels)
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

#' Wrapper class for server Perason correlation data
#'
#' @param data pearson r correlation coefficients
#' @param proj_labels the labels of the projections (columns) in the data
#' @param sig_labels the labels of the signatures (rows) in the data
#' @return a ServerPCorr object
ServerPCorr <- function(data, proj_labels, sig_labels) {
        .Object <- new("ServerPCorr", data=data, proj_labels=proj_labels,
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

    coord <- as.data.frame(t(p))

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
#' @param sigpm SigProjMatrix
#' @param sigs Signatures to subset form sigpm
#' @return Subsetted sigProjMatirx converted to JSON
sigProjMatrixToJSON <- function(sigpm, sigs) {

    sigpm <- sigpm[sigs,, drop=FALSE]
    sSPM <- ServerSigProjMatrix(unname(sigpm), colnames(sigpm), sigs)

    json <- toJSON(sSPM, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

    return(json)
}

#' convert perason correlation coeffcients between PCs and sgnatures into a JSON object
#' @param pc the pearson correlations matrx
#' @param sigs the signatures of interest
#' @return Subsetted pearson correlations converted to JSON
pearsonCorrToJSON <- function(pc, sigs) {

    pc <- pc[sigs,,drop=FALSE]
    cn <- c()
    for (i in 1:ncol(pc)) { cn <- c(cn, paste("PC", i)) }
    sPC <- ServerPCorr(unname(pc), cn, sigs)

    json <- toJSON(sPC, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

    return(json)

}

#' Converts the -log10(pvalues) of the consistency scores into a JSON object
#' @importFrom jsonlite toJSON
#' @param sigpmp SigProjMatrix p values
#' @param sigs Signatrues to subset from sigpmp
#' @return Subsetted sigProjMatrix_P converted to JSON
sigProjMatrixPToJSON <- function(sigpmp, sigs) {

    sigpmp <- as.matrix(sigpmp[sigs,, drop=FALSE])
    sPM <- ServerPMatrix(unname(sigpmp), colnames(sigpmp), sigs)

    json <- toJSON(sPM, force=TRUE, pretty=TRUE, auto_unbox=TRUE)

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

    url <- paste0("http://", host, ":", port, "/Results.html")
    message(paste("Navigate to", url, "in a browser to interact with the app."))

    if(browser){
        browseURL(url)
    }

    # Launch the server
    jug() %>%
      get("/Signature/Scores/(?<sig_name1>.*)", function(req, res, err) {
        sigMatrix <- object@sigMatrix
        name <- URLdecode(req$params$sig_name1)
        out <- "Signature does not exist!"
        if (name %in% colnames(sigMatrix)) {
          out <- FastProjectR:::sigScoresToJSON(sigMatrix[name])
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

        cls <- object@filterModuleData@ProjectionData@sigClusters
        cls <- cls$Computed

        out <- toJSON(cls, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/SigClusters/Meta", function(req, res, err) {

        cls <- object@filterModuleData@ProjectionData@sigClusters
        cls <- cls$Meta

        out <- toJSON(cls, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/(?<proj_name1>.*)/coordinates", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name1)
        out <- FastProjectR:::coordinatesToJSON(object@filterModuleData@ProjectionData@projections[[proj]]@pData)
        return(out)
      }) %>%
      get("/FilterGroup/projections/list", function(req, res, err) {
        proj_names <- names(object@filterModuleData@ProjectionData@projections)
        out <- toJSON(proj_names, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/SigProjMatrix/Normal", function(req, res, err) {

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isMeta, TRUE)
        sigs <- keys[!vals]
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@ProjectionData@sigProjMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/SigProjMatrix/Meta", function(req, res, err) {

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isMeta, TRUE)
        sigs <- keys[!vals]
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@ProjectionData@sigProjMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/SigProjMatrix_P/Normal", function(req, res, err) {

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isMeta, TRUE)
        sigs <- keys[!vals]
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@ProjectionData@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/SigProjMatrix_Pemp/Normal", function(req, res, err) {

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isMeta, TRUE)
        sigs <- keys[!vals]
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@ProjectionData@emp_pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/SigProjMatrix_P/Meta", function(req, res, err) {

        sigs <- colnames(object@metaData)

        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@ProjectionData@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/Tree/SigProjMatrix_P/Normal", function(req, res, err) {

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isMeta, TRUE)
        sigs <- keys[!vals]
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@TreeProjectionData@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/Tree/SigProjMatrix_P/Meta", function(req, res, err) {

        sigs <- colnames(object@metaData)
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleData@TreeProjectionData@pMatrix, sigs)

        return(out)
      }) %>%
      get("/FilterGroup/(?<proj_name2>.*)/clusters/(?<cluster_procedure>.*)/(?<param>.*)", function(req, res, err) {
        # projData <- object@projData

        proj <- URLdecode(req$params$proj_name2)
        method <- URLdecode(req$params$cluster_procedure)
        param <- as.numeric(URLdecode(req$params$param))

        clust <- FastProjectR:::cluster(object@filterModuleData@ProjectionData@projections[[proj]], method, param)
        out <- FastProjectR:::clusterToJSON(clust)
        return(out)
      }) %>%
      get("/FilterGroup/Tree/List", function(req, res, err) {

        ## all Trees have the same adjacency matrix, so we can use the first one
        W <- object@filterModuleData@TreeProjectionData@projections[[1]]@adjMat

        return(toJSON(W))
      }) %>%
      get("/FilterGroup/(?<proj_name3>.*)/Tree/Points", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name3)

        C <- object@filterModuleData@TreeProjectionData@projections[[proj]]@vData

        return(toJSON(C))
      }) %>%
      get("/FilterGroup/(?<proj_name4>.*)/Tree/Projection", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name4)

        out <- FastProjectR:::coordinatesToJSON(object@filterModuleData@TreeProjectionData@projections[[proj]]@pData)

        return(out)
      }) %>%
      get("/FilterGroup/list", function(req, res, err) {
        filters <- vapply(object@filterModuleList, function(x) {
          return(x@filter)
        }, "")
        return(toJSON(filters))
      }) %>%
      get("/FilterGroup/(?<pc_num1>.*)/Loadings/Positive", function(req, res, err) {
        pcnum <- as.numeric(URLdecode(req$params$pc_num1))

        l <- object@filterModuleData@PCAnnotatorData@loadings[,pcnum]

        posl <- l[l >= 0]
        sumposl <- sum(posl)
        nposl <- vapply(posl, function(x) x / sumposl, 1.0)

        nposl <- sort(nposl, decreasing = TRUE)

        js1 <- toJSON(with(stack(nposl), tapply(values, ind, c, simplify = FALSE)))

          return(js1)

      }) %>%
      get("/FilterGroup/(?<pc_num6>.*)/Loadings/Negative", function(req, res, err) {
        pcnum <- as.numeric(URLdecode(req$params$pc_num6))

        l <- object@filterModuleData@PCAnnotatorData@loadings[,pcnum]

        negl <- l[l < 0]
        sumneg1 <- sum(negl)
        nnegl <- vapply(negl, function(x) x / sumneg1, 1.0)

        nnegl <- sort(nnegl, decreasing = TRUE)

        js2 <- toJSON(with(stack(nnegl),
                           tapply(values, ind, c, simplify = FALSE)))

        return(js2)

      }) %>%
      get("/FilterGroup/PCA/Coordinates", function(req, res, err) {

        pc <- object@filterModuleData@PCAnnotatorData@fullPCA
        out <- FastProjectR:::coordinatesToJSON(pc)

        return(out)
      }) %>%
      get("/FilterGroup/PCVersus/(?<pc_num3>.*)/(?<pc_num4>.*)", function(req, res, err) {

        pc1 <- as.numeric(URLdecode(req$params$pc_num3))
        pc2 <- as.numeric(URLdecode(req$params$pc_num4))

        pcdata1 <- object@filterModuleData@PCAnnotatorData@fullPCA[pc1, ]
        pcdata2 <- object@filterModuleData@PCAnnotatorData@fullPCA[pc2, ]


        ret <- cbind(pcdata1, pcdata2)
        coord <- apply(unname(ret), 1, as.list)
        names(coord) <- rownames(ret)

        return(toJSON(coord, force = TRUE, auto_unbox = TRUE))
      }) %>%
      get("/FilterGroup/PearsonCorr/Normal", function(req, res, err) {

          signatures <- object@sigData
          keys <- vapply(signatures, function(x) x@name, "")
          vals <- vapply(signatures, function(x) x@isMeta, TRUE)
          sigs <- keys[!vals]
          pc <- object@filterModuleData@PCAnnotatorData@pearsonCorr

        return(FastProjectR:::pearsonCorrToJSON(pc, sigs))
      }) %>%
      get("/FilterGroup/PearsonCorr/Meta", function(req, res, err) {

          sigs <- colnames(object@metaData)
          numericMeta <- vapply(sigs,
                                function(x) is.numeric(object@metaData[[x]]),
                                FUN.VALUE = TRUE)
          sigs <- sigs[numericMeta]
          pc <- object@filterModuleData@PCAnnotatorData@pearsonCorr

        return(FastProjectR:::pearsonCorrToJSON(pc, sigs))
      }) %>%
      get("/FilterGroup/PearsonCorr/list", function(req, res, err) {

          pc <- object@filterModuleData@PCAnnotatorData@pearsonCorr
          pcnames <- seq(ncol(pc))
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

        data <- object@exprData
        gene_name <- URLdecode(req$params$gene_name2)

        result <- toJSON(
                         as.list(data[gene_name,]),
                         force=TRUE, pretty=TRUE, auto_unbox=TRUE
                         )

        return(result)

      }) %>%
      get("/SessionInfo", function(req, res, err) {

        info <- list()

        if (.hasSlot(object, "name") && !is.null(object@name)) {
            info["name"] <- object@name
        } else {
            info["name"] <- ""
        }

        W <- object@filterModuleData@TreeProjectionData
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
      serve_static_files(
                         root_path=system.file("html_output", package="FastProjectR")
                        ) %>%
      simple_error_handler_json() %>%
      serve_it(host=host, port=port)

}
