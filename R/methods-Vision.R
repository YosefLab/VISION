#' Initializes a new VISION object.
#'
#' @import logging
#'
#' @param data expression data - can be one of these: \itemize{
#' \item numeric matrix
#' \item ExpressionSet object
#' \item SummzrizedExperiment object (or extending classes)
#' }
#' @param signatures list of file paths to signature files (.gmt or .txt) or a
#' list of Signature objects
#' @param housekeeping vector of gene names
#' @param meta data file with cell meta data (.txt), or a
#' data.frame with meta information. Note that rows should match samples in the
#' data, and columns should either be factors or numerics.
#' @param nomodel if TRUE, no fnr curve calculated and all weights equal to 1.
#' Else FNR and weights calculated. [Default:TRUE]
#' @param projection_genes name of method ('threshold' or 'fano') or list of
#' genes to use when computing projections.
#' @param min_signature_genes Minimum number of genes required to compute a
#' signature
#' @param weights Precomputed weights for each coordinate. Normally computed
#' from the FNR curve.
#' @param threshold Threshold to apply for the threshold filter. Number of cells or
#' proportion of cells (if less than 1)
#' @param perm_wPCA If TRUE, apply permutation WPCA to calculate significant
#' number of PCs. Else not. Default FALSE.
#' @param sig_norm_method Method to apply to normalize the expression matrix
#' before calculating signature scores. Valid options are:
#' "znorm_columns" (default), "none", "znorm_rows", "znorm_rows_then_columns",
#' or "rank_norm_columns"
#' @param sig_score_method Method to apply when calculating signature scores.
#' Either "naive" (default) or "weighted_avg"
#' @param pool indicates whether or not to create supercells. Acceptable values
#' are TRUE, FALSE, or 'auto', the last of which is the default and enables
#' pooling if there are more than 15000 cells.
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @param cluster_variable variable to use to denote clusters
#' @param latentSpace latent space for expression data. Numeric matrix or dataframe
#' with dimensions CELLS x COMPONENTS
#' @param latentTrajectory trajectory to model cell progression.  Wrapped result
#' of a trajectory inference by the dynverse/dynwrap library
#' @param metaData a dataframe storing all metadata (can be either factor or numerical),
#' where the rows are cells and columns are individual meta data items
#' @param projection_methods a character list of which projection methods to apply. Can be: \itemize{
#'    \item tSNE10 (tSNE with perplexity 10)
#'    \item tSNE30 (tSNE with perplexity 30)
#'    \item ICA
#'    \item ISOMap
#'    \item RBFPCA
#'}
#' By default will perform tSNE and PCA on the data.
#' @param name a name for the sample - shown on the output report
#' @return A VISION object
#' @rdname VISION-class
#' @export
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
#'
#' #create 20 signatures of 25 genes each
#' sigs <- lapply(1:20, function(i) {
#' sigData <- sign(rnorm(25))
#' names(sigData) <- paste0("gene",sample(1:100,25))
#' return(createGeneSignature(name = paste0("sig",i),
#'                                  sigData = sigData))
#' })
#'
#' fp <- vision(data = expMat,
#'                      signatures = sigs,
#'                      housekeeping = hkg)
setMethod("Vision", signature(data = "matrixORSparse"),
            function(data, signatures, housekeeping=NULL,
                    unnormalizedData = NULL, meta=NULL, nomodel=TRUE,
                    projection_genes=c("fano"), min_signature_genes=5,
                    weights=NULL, threshold=.05, perm_wPCA=FALSE,
                    projection_methods = c("tSNE30"),
                    sig_norm_method = c("znorm_columns", "none", "znorm_rows",
                                        "znorm_rows_then_columns",
                                        "rank_norm_columns"),
                    sig_score_method=c("naive", "weighted_avg"),
                    pool="auto", cellsPerPartition=100, name=NULL,
                    cluster_variable = "", latentSpace = NULL,
                    latentTrajectory = NULL) {

            .Object <- new("Vision")

            rownames(data) <- toupper(rownames(data))
            .Object@initialExprData <- data
            .Object@exprData <- data

            if (!is.null(unnormalizedData)){

                if (is.data.frame(unnormalizedData)){
                    unnormalizedData <- data.matrix(unnormalizedData)
                }
                # unnormalizedData might have more genes than exprData
                # and it might have more cells than exprData
                rownames(unnormalizedData) <- toupper(rownames(unnormalizedData))
                HAS_CORRECT_CELLS <- length(setdiff(
                                           colnames(.Object@exprData),
                                           colnames(unnormalizedData)
                                           )) == 0
                if (!HAS_CORRECT_CELLS) {
                    stop("unnormalizedData must have a column for all cells in data. colnames(unnormalizedData) must contain all labels in colnames(data)")
                }

                HAS_CORRECT_GENES <- length(setdiff(
                                           rownames(.Object@exprData),
                                           rownames(unnormalizedData)
                                           )) == 0
                if (!HAS_CORRECT_GENES) {
                    stop("unnormalizedData must have a row for all genes in data. rownames(unnormalizedData) must contain all labels in rownames(data)")
                }

                if (any(unnormalizedData < 0)) {
                    stop("Negative values in unnormalizedData. unnormalizedData should be either counts or scaled counts and should therefore have no negative values")
                }

                .Object@unnormalizedData <- unnormalizedData
                .Object@initialUnnormalizedData <- unnormalizedData

            } else {
                .Object@unnormalizedData <- .Object@exprData
                .Object@initialUnnormalizedData <- .Object@initialExprData
            }

            if (is.null(housekeeping)) {
                .Object@housekeepingData <- character()
                .Object@nomodel = TRUE
            } else {
                .Object@housekeepingData <- vapply(housekeeping,
                                                   toupper, "",
                                                   USE.NAMES = FALSE)
            }

            if (is.list(signatures)) {
                .Object@sigData <- signatures
                names(.Object@sigData) <- vapply(.Object@sigData,
                                            function(x){x@name}, "")
            } else if (is.character(signatures)) {
                .Object@sigData <- readSignaturesInput(signatures)
            } else {
                stop("signatures must be paths to signature files or list of
                    Signature objects")
            }

            .Object@sigData <- processSignatures(.Object@sigData, rownames(.Object@exprData), min_signature_genes)

            if (!is.null(meta)) {
                if(is.matrix(meta)){
                    meta <- as.data.frame(meta)
                }
                if(is.data.frame(meta)) {
                    sampleLabels <- colnames(.Object@exprData)

                    common <- intersect(row.names(meta), sampleLabels)
                    if (length(common) != length(sampleLabels)){
                        stop("Provided meta data dataframe must have same sample labels as the expression matrix")
                    }

                    # Convert strings to factors if less than 20 unique
                    metaVars <- colnames(meta)
                    for(var in metaVars){
                        vals <- meta[, var]
                        if (is.character(vals)){
                            n_unique <- length(unique(vals))
                            if (n_unique <= 20){
                                meta[, var] <- as.factor(vals)
                            } else {
                                meta[, var] <- NULL
                                message(paste0("Dropping '", var, "' from meta data as it is of type 'character' and has more than 20 unique values.  If you want to include this meta data variable, convert it to a factor before providing the data frame to Vision"))
                            }
                        }

                        if (is.logical(vals)){
                            meta[, var] <- as.factor(vals)
                        }
                    }

                    .Object@metaData <- meta[sampleLabels, , drop = FALSE]
                    .Object@initialMetaData <- .Object@metaData
                } else {
                    stop("meta input argument should be a matrix or dataframe")
                }
            } else {
                .Object@metaData <- data.frame(
                                        row.names = colnames(.Object@exprData)
                                    )
                .Object@initialMetaData <- .Object@metaData
            }

            if (!is.null(weights)) {
                .Object@weights <- weights
            }

            if (!.Object@nomodel) {
                .Object@nomodel <- nomodel
            }
            .Object@projection_genes <- vapply(projection_genes,
                                               toupper, "",
                                               USE.NAMES = FALSE)

            if (threshold < 1) {
                num_samples <- ncol(.Object@exprData)
                threshold <- round(threshold * num_samples)
            }

            .Object@threshold <- threshold
            .Object@sig_norm_method <- match.arg(sig_norm_method)
            .Object@sig_score_method <- match.arg(sig_score_method)
            .Object@perm_wPCA <- perm_wPCA

            valid_projections <- c("tSNE10", "tSNE30", "ICA", "ISOMap", "RBFPCA")
            check <- sapply(projection_methods, function(x) x %in% valid_projections)
            if (! all(check)) {
                stop("Bad value in 'projection_methods'. Please choose from tSNE10, tSNE30, ICA, ISOMap, or RBFPCA.")
            }
            .Object@projection_methods <- projection_methods

            LOTS_OF_CELLS <- ncol(.Object@exprData) > 15000

            if (is.character(pool))
            {
                pool <- tolower(pool)
                if (pool == 'auto')
                {
                    if(LOTS_OF_CELLS) {
                        pool <- TRUE
                        message(paste(
                          "Over 15000 input cells detected.  Enabling micropooling with max",
                          cellsPerPartition, "cells per pool."
                          )
                        )
                    } else {
                        pool <- FALSE
                    }
                } else {
                    stop("Bad value for 'pool' argument")
                }
            }

            if (LOTS_OF_CELLS && !pool) {
                message(paste(
                      "Warning: Input data consists of",
                      ncol(.Object@exprData),
                      "cells and pool=FALSE.  It is recommend to set pool=TRUE when running on large numbers of cells"
                      )
                )
            }

            .Object@pool = pool
            .Object@cellsPerPartition = cellsPerPartition

            if (!is.null(name)) {
                .Object@name <- name
            }

            if (!is.null(latentSpace)) {
                if (is.data.frame(latentSpace)){
                    latentSpace <- data.matrix(latentSpace)
                }

                sample_names <- colnames(.Object@exprData)
                common <- intersect(sample_names, rownames(latentSpace))

                if (length(common) != nrow(latentSpace)){
                    stop("Supplied coordinates for latentSpace must have rowlabels that match sample/cell names")
                }

                latentSpace <- latentSpace[colnames(.Object@exprData), ]
                colnames(latentSpace) <- NULL

                .Object@latentSpace <- latentSpace
                .Object@initialLatentSpace <- latentSpace
            }

            if (!is.null(latentTrajectory)) {

                if (!(
                      "milestone_network" %in% names(latentTrajectory) &&
                      "progressions" %in% names(latentTrajectory)
                  )){
                    stop("latentTrajectory must be a wrapped method using the dynverse/dynmethods library.")
                }

                .Object@latentTrajectory <- Trajectory(latentTrajectory)

                sample_names <- colnames(.Object@exprData)
                if (length(
                        setdiff(
                            sample_names,
                            rownames(.Object@latentTrajectory@progressions)
                            )
                        ) > 0) {
                    stop("Supplied progressions for latentTrajectory must have cell_ids that match sample/cell names")
                }

                newMeta <- createTrajectoryMetaData(.Object@latentTrajectory)
                newMeta <- newMeta[rownames(.Object@metaData), ]

                .Object@metaData <- cbind(.Object@metaData, newMeta)
                .Object@initialMetaData <- .Object@metaData
            }

            .Object@cluster_variable <- cluster_variable

            return(.Object)
    }
)

#' @rdname VISION-class
#' @export
setMethod("Vision", signature(data = "data.frame"),
            function(data, ...) {
                data <- data.matrix(data)
            return(Vision(data, ...))
            }
)

#' @rdname VISION-class
#' @export
setMethod("Vision", signature(data = "ExpressionSet"),
            function(data, ...) {
              if (!requireNamespace("Biobase", quietly = TRUE)){
                  stop("Package \"Biobase\" needed to load this data object.  Please install it.",
                       call. = FALSE)
              }
            return(Vision(Biobase::exprs(data), ...))
            }
)

#' @rdname VISION-class
#' @export
setMethod("Vision", signature(data = "SummarizedExperiment"),
          function(data, ...) {

              if (!requireNamespace("SummarizedExperiment", quietly = TRUE)){
                  stop("Package \"SummarizedExperiment\" needed to load this data object.  Please install it.",
                       call. = FALSE)
              }
            return(Vision(SummarizedExperiment::assay(data), ...))
          }
)

#' Main entry point for running VISION Analysis
#'
#' The main analysis function. Runs the entire VISION analysis pipeline
#' and returns a VISION object populated with the result,
#'
#' @export
#' @aliases analyze
#' @param object VISION object
#' @return VISION object
#'
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
#'
#' #create 20 signatures of 25 genes each
#' sigs <- lapply(1:20, function(i) {
#' sigData <- sign(rnorm(25))
#' names(sigData) <- paste0("gene",sample(1:100,25))
#' return(createGeneSignature(name = paste0("sig",i),
#'                                  sigData = sigData))
#' })
#'
#' fp <- VISION(data = expMat,
#'                      housekeeping = hkg,
#'                      signatures = sigs)
#'
#' ## analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- analyze(fp)
#' }
setMethod("analyze", signature(object="Vision"),
            function(object) {
    message("Beginning Analysis")

    if (object@cluster_variable == "") {
        object <- clusterCells(object)
    }

    if (object@pool) {
        object <- poolCells(object)
    }

    object <- filterData(object)
    object <- convertToDense(object)

    object <- calcWeights(object)

    # Populates @sigScores
    object <- calcSignatureScores(object)

    # Populates @latentSpace
    if (all(dim(object@latentSpace) == c(1, 1))) {
        object <- computeLatentSpace(object)
    }

    # Populates @Projections
    object <- generateProjections(object)

    # Populates @TrajectoryProjections
    if (!is.null(object@latentTrajectory)) {

        object@TrajectoryProjections <- generateTrajectoryProjections(
                                            object@latentTrajectory
                                        )
    }

    message("Computing background distribution for signature scores...")
    signatureBackground <- calculateSignatureBackground(object, num = 3000)

    # Populates @SigConsistencyScores
    object <- analyzeSpatialCorrelations(object, signatureBackground)

    # Populates @TrajectoryConsistencyScores
    if (!is.null(object@latentTrajectory)) {
        object <- analyzeTrajectoryCorrelations(object, signatureBackground)
    }

    # Populates @ClusterSigScores
    object <- clusterSigScores(object)

	# Add Numerical MetaData to @sigScores
    object <- addNumericalMetaData(object)

    # Populates #PCAnnotatorData
    object <- calculatePearsonCorr(object)


    message("Analysis Complete!")

    return(object)
})

#' Add a set of projection coordinates
#'
#' @export
#' @param object VISION object
#' @param name Name of the projection
#' @param coordinates numeric matrix or data.frame. Coordinates of each
#' sample in the projection (NUM_SAMPLES x NUM_COMPONENTS)
#' @return VISION object
setMethod("addProjection", signature(object = "Vision"),
            function(object, name, coordinates) {

    if (is(coordinates, "data.frame")){
        coordinates <- as.matrix(coordinates)
    }

    # Verify that projection coordinates are correct
    samples <- object@exprData

    SAME_SIZE <- ncol(samples) == nrow(coordinates)
    SAME_NAMES <- setequal(colnames(samples), rownames(coordinates))

    if (!SAME_SIZE || !SAME_NAMES){
        stop("Supplied coordinates must have rowlabels that match sample/cell names")
    }

    if (dim(coordinates)[2] != 2){
        stop("Projection must have exactly 2 components")
    }

    object@inputProjections[[name]] <- coordinates

    return(object)
})

#' Save the VISION object as an .RDS file and view the results on a
#' localhost
#'
#' Save the results object as an RDS file for future use, and launch a local
#' server to explore the results with a browser.
#'
#' @param fpout VISION object
#' @param ofile the path to save the object in. If NULL, the object is saved
#' in the working directory [default:NULL]
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @return the path of the saved file
#' @aliases saveAndViewResults
#' @export
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
#'
#' #create 20 signatures of 25 genes each
#' sigs <- lapply(1:20, function(i) {
#' sigData <- sign(rnorm(25))
#' names(sigData) <- paste0("gene",sample(1:100,25))
#' return(createGeneSignature(name = paste0("sig",i),
#'                                  sigData = sigData))
#' })
#'
#' fp <- VISION(data = expMat,
#'                      housekeeping = hkg,
#'                      signatures = sigs)
#'
#' ## analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- analyze(fp)
#' saveAndViewResults(fp.out)
#' }
setMethod("saveAndViewResults", signature(fpout="Vision"),
          function(fpout, ofile=NULL, port=NULL, host=NULL, browser=TRUE, name=NULL) {
            if(is.null(ofile)) {
              i <- 1
              ofile <- paste0("./fpout", i, ".rds")
              while (file.exists(ofile)) {
                i <- i+1
                ofile <- paste0("./fpout", i, ".rds")
              }
            }

            saveRDS(fpout, file=ofile)
            viewResults(fpout, port, host, browser, name)
            return(ofile)
          })

#' View results of analysis
#'
#' launch a local server to explore the results with a browser.
#'
#' @param object VISION object or path to a file containing such an
#' object (saved using saveAndViewResults, or directly using saveRDS)
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @param browser Whether or not to launch the browser automatically (default=TRUE)
#' @param name Name for the sample - is shown at the top of the output report
#' @aliases viewResults
#' @return None
#' @export
#' @rdname viewResults
#' @examples
#' expMat <- matrix(rnorm(200000), nrow=500)
#' rownames(expMat) <- paste0("gene",1:500)
#'
#' # choose housekeeping genes
#' hkg <- housekeeping$default
#'
#' #create 20 signatures of 25 genes each
#' sigs <- lapply(1:20, function(i) {
#' sigData <- sign(rnorm(25))
#' names(sigData) <- paste0("gene",sample(1:100,25))
#' return(createGeneSignature(name = paste0("sig",i),
#'                                  sigData = sigData))
#' })
#'
#' fp <- VISION(data = expMat,
#'                      housekeeping = hkg,
#'                      signatures = sigs)
#'
#' ## analyze requires actual non-random data to run properly
#' \dontrun{
#' fp.out <- analyze(fp)
#' viewResults(fp.out)
#' }
setMethod("viewResults", signature(object="Vision"),
          function(object, port=NULL, host=NULL, browser=TRUE, name=NULL) {

            if (!is.null(name)) {
                object@name <- name
            }

            versionCheck(object)

            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            launchServer(object, port, host, browser)
          })

#' @rdname viewResults
#' @export
setMethod("viewResults", signature(object="character"),
          function(object, port=NULL, host=NULL, browser=TRUE, name=NULL) {
            fpo <- readRDS(object)
            if(!is(fpo, "Vision")){
              stop("loaded object not a valid Vision object")
            }
            viewResults(fpo, port, host, browser, name)
          })

#' create new VISION object from a subset of the data in an existing one
#' @param fp the VISION object to subset
#' @param subset the indices of the samples to keep
#' @return a new VISION object with the new data and the same analysis
#' parameters
createNewFP <- function(fp, subset) {
    .Object <- new("Vision")
    nexpr <- fp@initialExprData[,subset]
    rownames(nexpr) <- toupper(rownames(nexpr))
    .Object@initialExprData <- nexpr
    .Object@exprData <- nexpr

    .Object@housekeepingData <- fp@housekeepingData
    .Object@sigData <- fp@sigData

    .Object@metaData <- lapply(fp@metaData, function(sigscore) {
        sigscore@scores <- sigscore@scores[subset]
        return(sigscore)
    })

    if (!all(dim(fp@weights) == c(1, 1))){
        .Object@weights <- fp@weights[, subset]
    }
    .Object@projection_genes <- fp@projection_genes
    .Object@threshold <- fp@threshold
    .Object@sig_norm_method <- fp@sig_norm_method
    .Object@sig_score_method <- fp@sig_score_method
    .Object@perm_wPCA = fp@perm_wPCA
    .Object@pool = fp@pool
    .Object@cellsPerPartition = fp@cellsPerPartition

    return(.Object)
}
