#' Wrapper class for ExpressionData object for JSON.
#'
#' @param data Expression data
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

#' Converts list of signature ranks to JSON
#' @importFrom jsonlite toJSON
#' @param ss single-column dataframe with ranks for a single signature
#' @return Signature ranks as JSON, with names of each entry that of list names
sigRanksToJSON <- function(ss) {

    s <- as.list(rank(ss[[1]]))
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
#' @param object FastProject object or path to a file containing such an
#' object (saved using saveAndViewResults, or directly using saveRDS)
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
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
      get("/Signature/ListPrecomputed", function(req, res, err){
        signatures <- object@sigData
        keys <- lapply(signatures, function(x) x@name)
        vals <- lapply(signatures, function(x) x@isPrecomputed)
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
      get("/Signature/Ranks/(?<sig_name3>.*)", function(req, res, err) {
        sigMatrix <- object@sigMatrix
        name <- URLdecode(req$params$sig_name3)
        out <- "Signature does not exist!"
        if (name %in% colnames(sigMatrix)) {
          out <- FastProjectR:::sigRanksToJSON(sigMatrix[name])
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
            if(sig@isPrecomputed) {
                stop("Can't get expression for precomputed signature")
            }
            genes = names(sig@sigDict)
            expMat = object@exprData@data
            return(FastProjectR:::expressionToJSON(expMat, genes, zscore=TRUE))
        }
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_group1>.*)/SigClusters/Normal", function(req, res, err) {
        filter <- URLdecode(req$params$filter_group1)
        if (filter == "1") {
          filter = 1
        }
        cls <- object@filterModuleList[[filter]]@ProjectionData@sigClusters
        # cls <- object@sigClusters[[filter]]
        cls <- cls$Computed

        out <- toJSON(cls, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_group2>.*)/SigClusters/Precomputed", function(req, res, err) {
        filter <- URLdecode(req$params$filter_group2)
        if (filter == "1") {
          filter = 1
        }
        cls <- object@filterModuleList[[filter]]@ProjectionData@sigClusters
        # cls <- object@sigClusters[[filter]]
        cls <- cls$Precomputed

        out <- toJSON(cls, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name1>.*)/(?<proj_name1>.*)/coordinates", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name1)
        proj <- URLdecode(req$params$proj_name1)
        out <- FastProjectR:::coordinatesToJSON(object@filterModuleList[[filter]]@ProjectionData@projections[[proj]]@pData)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name1p5>.*)/projections/list", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name1p5)
        proj_names <- names(object@filterModuleList[[filter]]@ProjectionData@projections)
        out <- toJSON(proj_names, auto_unbox=TRUE)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name2>.*)/SigProjMatrix/Normal", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name2)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
          sigs <- keys[!vals]
        out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@ProjectionData@sigProjMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name3>.*)/SigProjMatrix/Precomputed", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name3)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
        sigs <- keys[!vals]
          out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@ProjectionData@sigProjMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name4>.*)/SigProjMatrix_P/Normal", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name4)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
        sigs <- keys[!vals]
          out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@ProjectionData@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name22>.*)/SigProjMatrix_Pemp/Normal", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name22)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
        sigs <- keys[!vals]
          out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@ProjectionData@emp_pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name5>.*)/SigProjMatrix_P/Precomputed", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name5)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
        sigs <- keys[vals]

          out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@ProjectionData@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name20>.*)/Tree/SigProjMatrix_P/Normal", function(req, res, err) {
        # projData <- object@projData
        filter <- URLdecode(req$params$filter_name20)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
        sigs <- keys[!vals]
          out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@TreeProjectionData@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name21>.*)/Tree/SigProjMatrix_P/Precomputed", function(req, res, err) {
        # projData <- object@projData
        filter <- URLdecode(req$params$filter_name21)

        signatures <- object@sigData
        keys <- vapply(signatures, function(x) x@name, "")
        vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
        sigs <- keys[vals]
          out <- FastProjectR:::sigProjMatrixToJSON(object@filterModuleList[[filter]]@TreeProjectionData@pMatrix, sigs)
        # out <- sigProjMatrixPToJSON(projData[[filter]]@pMatrix, sigs)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name6>.*)/(?<proj_name2>.*)/clusters/(?<cluster_procedure>.*)/(?<param>.*)", function(req, res, err) {
        # projData <- object@projData

        filter <- URLdecode(req$params$filter_name6)
        proj <- URLdecode(req$params$proj_name2)
        method <- URLdecode(req$params$cluster_procedure)
        param <- as.numeric(URLdecode(req$params$param))

        clust <- FastProjectR:::cluster(object@filterModuleList[[filter]]@ProjectionData@projections[[proj]], method, param)
          # clust = cluster(projData[[filter]]@projections[[proj]], method, param)
        out <- FastProjectR:::clusterToJSON(clust)
        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name7>.*)/genes", function(req, res, err) {
        # projData <- object@projData
        filter <- URLdecode(req$params$filter_name7)

        out <- toJSON(object@filterModuleList[[filter]]@genes)
        # out <- toJSON(projData[[filter]]@genes)

        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name8>.*)/Tree/List", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name8)

        ## all Trees have the same adjacency matrix, so we can use the first one
        W <- object@filterModuleList[[filter]]@TreeProjectionData@projections[[1]]@adjMat

        return(toJSON(W))
      }) %>%
      get("/FilterGroup/(?<filter_name18>.*)/(?<proj_name3>.*)/Tree/Points", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name3)
        filter <- URLdecode(req$params$filter_name18)

        C <- object@filterModuleList[[filter]]@TreeProjectionData@projections[[proj]]@vData

        return(toJSON(C))
      }) %>%
      get("/FilterGroup/(?<filter_name19>.*)/(?<proj_name4>.*)/Tree/Projection", function(req, res, err) {
        proj <- URLdecode(req$params$proj_name4)
        filter <- URLdecode(req$params$filter_name19)

        out <- FastProjectR:::coordinatesToJSON(object@filterModuleList[[filter]]@TreeProjectionData@projections[[proj]]@pData)

        return(out)
      }) %>%
      get("/FilterGroup/list", function(req, res, err) {
        filters <- vapply(object@filterModuleList, function(x) {
          return(x@filter)
        }, "")
        return(toJSON(filters))
      }) %>%
      get("/FilterGroup/(?<filter_name10>.*)/(?<pc_num1>.*)/Loadings/Positive", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name10)
        pcnum <- as.numeric(URLdecode(req$params$pc_num1))

        l <- object@filterModuleList[[filter]]@PCAnnotatorData@loadings[,pcnum]

        posl <- l[l >= 0]
        sumposl <- sum(posl)
        nposl <- vapply(posl, function(x) x / sumposl, 1.0)

        nposl <- sort(nposl, decreasing=T)

        js1 <- toJSON(with(stack(nposl), tapply(values, ind, c, simplify=F)))

          return(js1)

      }) %>%
      get("/FilterGroup/(?<filter_name16>.*)/(?<pc_num6>.*)/Loadings/Negative", function(req, res, err) {
        filter <- URLdecode(req$params$filter_name16)
        pcnum <- as.numeric(URLdecode(req$params$pc_num6))

        l <- object@filterModuleList[[filter]]@PCAnnotatorData@loadings[,pcnum]

        negl <- l[l < 0]
        sumneg1 <- sum(negl)
        nnegl <- vapply(negl, function(x) x / sumneg1, 1.0)

        nnegl <- sort(nnegl, decreasing=T)

        js2 <- toJSON(with(stack(nnegl), tapply(values, ind, c, simplify=F)))

        return(js2)

      }) %>%
      get("/FilterGroup/(?<filter_name12>.*)/PCA/Coordinates", function(req, res, err) {
        # projData <- object@projData
        filter <- URLdecode(req$params$filter_name12)

        pc <- object@filterModuleList[[filter]]@PCAnnotatorData@fullPCA
        out <- FastProjectR:::coordinatesToJSON(pc)

        return(out)
      }) %>%
      get("/FilterGroup/(?<filter_name12>.*)/(?<sig_name5>.*)/(?<pc_num2>.*)/Coordinates", function(req, res, err) {
        # projData <- object@projData
        filter <- URLdecode(req$params$filter_name12)
        pcnum <- as.numeric(URLdecode(req$params$pc_num2))
        signame <- URLdecode(req$params$sig_name5)

        pc <- object@filterModuleList[[filter]]@PCAnnotatorData@fullPCA[pcnum,]
        # pc <- projData[[filter]]@fullPCA[pcnum,]
        ss <- object@sigMatrix[,signame]

        ret <- cbind(pc, ss)
        coord <- apply(unname(ret), 1, as.list)
        names(coord) <- rownames(ret)

        return(toJSON(coord, force=T, auto_unbox=T))
      }) %>%
      get("/FilterGroup/(?<filter_name13>.*)/PCVersus/(?<pc_num3>.*)/(?<pc_num4>.*)", function(req, res, err) {
        # projData <- object@projData
        filter <- URLdecode(req$params$filter_name13)
        pc1 <- as.numeric(URLdecode(req$params$pc_num3))
        pc2 <- as.numeric(URLdecode(req$params$pc_num4))

        pcdata1 <- object@filterModuleList[[filter]]@PCAnnotatorData@fullPCA[pc1,]
        # pcdata1 <- projData[[filter]]@fullPCA[pc1,]
        pcdata2 <- object@filterModuleList[[filter]]@PCAnnotatorData@fullPCA[pc2,]
        # pcdata2 <- projData[[filter]]@fullPCA[pc2,]


        ret <- cbind(pcdata1, pcdata2)
        coord <- apply(unname(ret), 1, as.list)
        names(coord) <- rownames(ret)

        return(toJSON(coord, force=T, auto_unbox=T))
      }) %>%
      get("/FilterGroup/(?<filter_name14>.*)/PearsonCorr/Normal", function(req, res, err) {
          filter <- URLdecode(req$params$filter_name14)

          signatures <- object@sigData
          keys <- vapply(signatures, function(x) x@name, "")
          vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
          sigs <- keys[!vals]
          pc <- object@filterModuleList[[filter]]@PCAnnotatorData@pearsonCorr

        return(FastProjectR:::pearsonCorrToJSON(pc, sigs))
      }) %>%
      get("/FilterGroup/(?<filter_name15>.*)/PearsonCorr/Precomputed", function(req, res, err) {
          filter <- URLdecode(req$params$filter_name15)

          signatures <- object@sigData
          keys <- vapply(signatures, function(x) x@name, "")
          vals <- vapply(signatures, function(x) x@isPrecomputed, TRUE)
          vals2 <- vapply(signatures, function(x) !x@isFactor, TRUE)
          sigs <- keys[vals & vals2]
          pc <- object@filterModuleList[[filter]]@PCAnnotatorData@pearsonCorr

        return(FastProjectR:::pearsonCorrToJSON(pc, sigs))
      }) %>%
      get("/FilterGroup/(?<filter_name17>.*)/PearsonCorr/list", function(req, res, err) {
          filter <- URLdecode(req$params$filter_name17)

          pc <- object@filterModuleList[[filter]]@PCAnnotatorData@pearsonCorr
          pcnames <- seq(ncol(pc))
          result <- toJSON(
                           pcnames,
                           force=TRUE, pretty=TRUE, auto_unbox=TRUE
                           )
          return(result)
      }) %>%
      get("/Expression/Genes/List", function(req, res, err) {

        data <- getExprData(object@exprData)
        genes = rownames(data)

        result <- toJSON(
                         genes,
                         force=TRUE, pretty=TRUE, auto_unbox=TRUE
                         )

        return(result)

      }) %>%
      get("/Expression/Gene/(?<gene_name2>.*)", function(req, res, err) {

        data <- getExprData(object@exprData)
        gene_name <- URLdecode(req$params$gene_name2)

        result <- toJSON(
                         as.list(data[gene_name,]),
                         force=TRUE, pretty=TRUE, auto_unbox=TRUE
                         )

        return(result)

      }) %>%
      post("/Analysis/Run/", function(req, res, err) {
        subset <- fromJSON(req$body)
        subset <- subset[!is.na(subset)]
        allData <- object@allData

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
