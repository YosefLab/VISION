#' Initializes a new VISION object.
#'
#' @import logging
#'
#' @param data expression data - can be one of these: \itemize{
#' \item numeric matrix or sparse matrix (GENES x CELLS)
#' \item data.frame (GENES x CELLS)
#' \item ExpressionSet object
#' \item SummzrizedExperiment object (or extending classes)
#' }
#' Expression data should be scaled and normalized, but not log-transformed.
#' @param signatures list of file paths to signature files (.gmt or .txt) or
#' Signature objects.  See the createGeneSignature(...) method for information
#' on creating Signature objects.
#' @param housekeeping vector of gene names
#' @param meta data.frame with meta-data for cells. Rows in this data.frame should correspond
#' with columns in the expression data matrix
#' @param nomodel if TRUE, no fnr curve calculated and all weights equal to 1.
#' Else FNR and weights calculated. [Default:TRUE]
#' @param projection_genes name of filtering method ('threshold' or 'fano') or list of
#' genes to use when computing projections.
#' @param min_signature_genes Signature that match less than this number of genes in the
#' supplied expression matrix are removed.
#' @param weights Precomputed weights for each coordinate. Normally computed
#' from the FNR curve.
#' @param threshold Threshold to apply when using the 'threshold' projection genes filter.
#' If greater than 1, this specifies the number of cells in which a gene must be detected
#' for it to be used when computing PCA. If less than 1, this instead specifies the proportion of cells needed
#' @param perm_wPCA If TRUE, apply permutation procedure to calculate significant
#' number of PCs when running PCA.  If FALSE (default), retain the top 30.
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
#' @param latentSpace latent space for expression data. Numeric matrix or dataframe
#' with dimensions CELLS x COMPONENTS
#' @param latentTrajectory trajectory to model cell progression.  Wrapped result
#' of a trajectory inference method by the dynverse/dynwrap library
#' @param projection_methods a character vector specifying which projection methods to apply. Can be: \itemize{
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
#' \dontrun{
#' expMat <- read.csv("expressionMatrix.csv", row.names=1)
#' meta <- read.csv("metaData.csv", row.names=1)
#'
#' sigs <- c("/path/to/signatures/msigdb_Hallmark.gmt",
#'           "/path/to/signatures/Enrichr/ChEA_2015.txt"
#'          )
#'
#'
#' vis <- Vision(data = expMat,
#'               signatures = sigs,
#'               meta = meta)
#' }
setMethod("Vision", signature(data = "matrixORSparse"),
            function(data, signatures=list(), housekeeping=NULL,
                    unnormalizedData = NULL, meta=NULL, nomodel=TRUE,
                    projection_genes=c("fano"), min_signature_genes=5,
                    weights=NULL, threshold=.05, perm_wPCA=FALSE,
                    projection_methods = c("tSNE30"),
                    sig_norm_method = c("znorm_columns", "none", "znorm_rows",
                                        "znorm_rows_then_columns",
                                        "rank_norm_columns"),
                    sig_score_method=c("naive", "weighted_avg"),
                    pool="auto", cellsPerPartition=100, name=NULL,
                    latentSpace = NULL, latentTrajectory = NULL, pools=list()) {

            .Object <- new("Vision")


            if (is.null(rownames(data))) {
                stop("rownames(data) = NULL. Expression matrix must have gene names as the rownames")
            }

            if (is.null(colnames(data))) {
                colnames(data) <- paste0("cell", seq(ncol(data)))
            }

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
                sigs <- lapply(signatures, function(sig){
                           if (is(sig, "Signature")){
                               return(sig)
                           } else {
                               return(readSignaturesInput(sig))
                           }
                })
                if (length(sigs) > 0){
                    .Object@sigData <- do.call(c, sigs)
                } else {
                    .Object@sigData <- sigs
                }
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
            if(length(projection_methods) == 0){
                projection_methods <- character()
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

            .Object@pools <- pools

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
#' and returns a VISION object populated with the result.
#'
#' @export
#' @aliases analyze
#' @param object VISION object
#' @return VISION object
#'
#' @examples
#' \dontrun{
#'
#' vis <- Vision(data = expMat, signatures = sigs)
#'
#' options(mc.cores=10) # Use 10 cores
#' vis <- analyze(vis)
#'
#' }
setMethod("analyze", signature(object="Vision"),
            function(object) {
    message("Beginning Analysis")

    if (object@cluster_variable == "") {
        object <- clusterCells(object)
    }

    if (object@pool || length(object@pools) > 0) {
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
    object <- analyzeLocalCorrelations(object, signatureBackground)

    # Populates @TrajectoryConsistencyScores
    if (!is.null(object@latentTrajectory)) {
        object <- analyzeTrajectoryCorrelations(object, signatureBackground)
    }

    # Populates @ClusterSigScores
    object <- clusterSigScores(object)

    # Populates #PCAnnotatorData
    object <- calculatePearsonCorr(object)


    message("Analysis Complete!")

    return(object)
})

#' Add a set of projection coordinates to use for visualization
#'
#' By default VISION will run tSNE on the latent space and use this
#' to display the cells in the output report.  However, with this
#' method, you may add additional two-dimensional coordinates for
#' inclusion in the output report. This is useful if you have previously
#' run tSNE (or any other visualization method) and wish to integrate
#' the results into VISION.
#'
#' @export
#' @aliases addProjection
#' @param object VISION object
#' @param name Name of the projection
#' @param coordinates numeric matrix or data.frame. Coordinates of each
#' sample in the projection (NUM_SAMPLES x NUM_COMPONENTS)
#' @return VISION object
#' @examples
#' \dontrun{
#'
#' # First create the VISION object
#' vis <- Vision(data = expMat, signatures = sigs)
#'
#' # Load and add an additional visualization
#' my_umap <- read.csv("umap_results.csv")
#' vis <- addProjection(vis, "UMAP", my_umap)
#'
#' # Run analysis
#' vis <- analyze(vis)
#'
#' # View results
#' viewResults(vis)
#'
#' }
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
#' This is just a convience function wrapping two function calls
#'
#'     saveAndViewResults(vis, 'vision_results.rds')
#'
#' is equivalent to:
#'
#'     saveRDS(vis, 'vision_results.rds')
#'     viewResults(vis)
#'
#' @param object VISION object
#' @param ofile the path to save the object in. If NULL, the object is saved
#' in the working directory [default:NULL]
#' @param port The port on which to serve the output viewer.  If omitted, a
#' random port between 8000 and 9999 is chosen.
#' @param host The host used to serve the output viewer. If omitted, "127.0.0.1"
#' is used.
#' @param browser Whether or not to launch the browser automatically (default=TRUE)
#' @param name Name for the sample - is shown at the top of the output report
#' @return the path of the saved file
#' @aliases saveAndViewResults
#' @export
#' @examples
#' \dontrun{
#'
#' vis <- Vision(data = expMat, signatures = sigs)
#'
#' options(mc.cores=10) # Use 10 cores
#' vis <- analyze(vis)
#'
#' saveAndViewResults(vis, 'vision_results.rds') # Saves and launches dynamic output report
#' }
setMethod("saveAndViewResults", signature(object = "Vision"),
          function(object, ofile=NULL, port=NULL, host=NULL,
                   browser=TRUE, name=NULL) {
            if (is.null(ofile)) {
              i <- 1
              ofile <- paste0("./vis", i, ".rds")
              while (file.exists(ofile)) {
                i <- i + 1
                ofile <- paste0("./vis", i, ".rds")
              }
            }

            saveRDS(object, file = ofile)
            viewResults(object, port, host, browser, name)
            return(ofile)
          })

#' View results of analysis
#'
#' launch a local server to explore the results with a browser.
#'
#' @param object VISION object or path to an RDS file containing such an
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
#' \dontrun{
#'
#' vis <- Vision(data = expMat, signatures = sigs)
#'
#' options(mc.cores=10) # Use 10 cores
#' vis <- analyze(vis)
#'
#' saveRDS(vis, 'vision_results.rds') # (Optional) Save results
#'
#' vis <- viewResults(vis) # Launches dynamic output report
#' }
setMethod("viewResults", signature(object = "Vision"),
          function(object, port=NULL, host=NULL, browser=TRUE, name=NULL) {

            if (!is.null(name)) {
                object@name <- name
            }

            versionCheck(object)

            if (length(object@sigData) == 0 && ncol(object@metaData) == 0) {
                stop("Error: This object contains no signature data.")
            }

            if (is.null(object@PCAnnotatorData)) {
                stop("Error: This object contains no PCAnnotatorData.")
            }

            if (!.hasSlot(object, "selections")) {
                object@selections <- list()
            }

            message("Launching the server...")
            message("Press exit or ctrl c to exit")
            object <- launchServer(object, port, host, browser)

            return(object)
          })

#' @rdname viewResults
#' @export
setMethod("viewResults", signature(object = "character"),
          function(object, port=NULL, host=NULL, browser=TRUE, name=NULL) {
            fpo <- readRDS(object)

            if (!is(fpo, "Vision")){
              stop("loaded object not a valid Vision object")
            }

            fpo <- viewResults(fpo, port, host, browser, name)
            return(fpo)
          })


#' Get saved selections
#'
#' Access saved groups of cell IDs defined while using the interactive output report
#'
#' This method allows you to retrieve saved selections later in R for downstream analyses
#'
#' Note:  In order for selections to correctly save when launching the report, the report
#'        must be run by storing the results back into the object.
#'
#' E.g.
#' \preformatted{vis <- viewResults(vis)}
#' and not
#' \preformatted{viewResults(vis)}
#'
#'
#' @param object VISION object
#' @return Named list of selections.  Each selection is a character vector of cell/pool IDs
#' @export
#' @aliases getSelections
#' @rdname getSelections
#' @examples
#' \dontrun{
#'
#' vis <- viewResults(vis)  # Selections saved while viewing results
#'
#' # Retrieve cell IDs for a selection group named 'interesting cells'
#' interestingCells <- getSelections(vis)[['interesting cells']]
#'
#' }
setMethod("getSelections", signature(object = "Vision"),
          function(object) {
              return(object@selections)
          })


#' Get 2D views of the expression data
#'
#' This method provides access to the 2d projections that are used
#' to display results in the output report
#'
#' @param object VISION object
#' @return List of matrix (Cells x 2)
#' @export
#' @aliases getProjections
#' @rdname getProjections
#' @examples
#' \dontrun{
#'
#' # After running 'analyze'
#' # Retrieve tSNE30 (tSNE with perplexity 30) and plot it
#'
#' tsne <- getProjections(vis)[["tSNE30"]]
#'
#' plot(tsne[, 1], tsne[, 2])
#'
#' # To see the names of available projections, just run:
#'
#' names(getProjections(vis))
#'
#' }
setMethod("getProjections", signature(object = "Vision"),
          function(object) {
              return(object@Projections)
          })


#' Get Latent Space
#'
#' Provides access to the latent space used for
#' local consistency analysis
#'
#' If a latent trajectory was supplied, access it by using \code{getLatentTrajectory}
#' instead
#'
#' @param object VISION object
#' @return the latent space as a matrix of dimension (Cells x Components)
#' @export
#' @aliases getLatentSpace
#' @rdname getLatentSpace
setMethod("getLatentSpace", signature(object = "Vision"),
          function(object) {
              return(object@latentSpace)
          })


#' Get Latent Trajectory
#'
#' Provides access to the latent trajectory used for
#' local consistency analysis
#'
#' If a latent space was supplied, access it by using \code{getLatentSpace}
#' instead
#'
#' @param object VISION object
#' @return Trajectory object
#' @export
#' @aliases getLatentTrajectory
#' @rdname getLatentTrajectory
#' @examples
#' \dontrun{
#'
#' trajectory <- getLatentTrajectory(vis)
#'
#' # MxM connectivity for network milestones
#' trajectory@adjMat
#'
#' # data.frame with the position of cells between milestones
#' # Columns are:
#' #    cell
#' #    from (milestone)
#' #    to (milestone)
#' #    position (0 to 1)
#' trajectory@progressions
#'
#' }
setMethod("getLatentTrajectory", signature(object = "Vision"),
          function(object) {
              return(object@latentTrajectory)
          })


#' Get Signature Scores
#'
#' Access to the signature scores computed by VISION
#'
#' @param object VISION object
#' @return Signature scores as a (Cells x Signature) matrix
#' @export
#' @aliases getSignatureScores
#' @rdname getSignatureScores
setMethod("getSignatureScores", signature(object = "Vision"),
          function(object) {
              return(object@sigScores)
          })


#' Get Signature Consistency Scores
#'
#' Access the local consistency scores computed for signatures
#'
#' Local consistency scores are calculated from the input latent
#' space (default's to PCA) or the input trajectory model (if provided)
#'
#' @param object VISION object
#' @return data.frame with columns 'Consistency', 'pValue', and 'FDR'
#' @export
#' @aliases getSignatureConsistency
#' @rdname getSignatureConsistency
setMethod("getSignatureConsistency", signature(object = "Vision"),
          function(object) {
              if (is.null(object@TrajectoryConsistencyScores)){
                  localData <- object@SigConsistencyScores
              } else {
                  localData <- object@TrajectoryConsistencyScores
              }

              # Remove meta-variables
              metaVars <- colnames(object@metaData)

              fdr <- localData@FDR[, 1]
              pVals <- localData@pValue[, 1]
              consistency <- localData@Consistency[, 1]

              sigs <- setdiff(names(pVals), metaVars)

              pVals <- pVals[sigs]
              consistency <- consistency[sigs]
              fdr <- fdr[sigs]

              out <- data.frame(
                  Consistency = consistency,
                  pValue = pVals,
                  FDR = fdr,
                  row.names = sigs
              )

              out <- out[order(out$pValue, out$Consistency * -1), ]

              return(out)

          })


#' Get MetaData Consistency Scores
#'
#' Access the local consistency scores computed for meta-data variables
#'
#' Local consistency scores are calculated from the input latent
#' space (default's to PCA) or the input trajectory model (if provided)
#'
#' @param object VISION object
#' @return data.frame with columns 'Consistency', 'pValue', and 'FDR'
#' @export
#' @aliases getMetaConsistency
#' @rdname getMetaConsistency
setMethod("getMetaConsistency", signature(object = "Vision"),
          function(object) {
              if (is.null(object@TrajectoryConsistencyScores)){
                  localData <- object@SigConsistencyScores
              } else {
                  localData <- object@TrajectoryConsistencyScores
              }

              # Remove non meta-variables
              metaVars <- colnames(object@metaData)

              fdr <- localData@FDR[, 1]
              pVals <- localData@pValue[, 1]
              consistency <- localData@Consistency[, 1]

              sigs <- intersect(names(pVals), metaVars)

              pVals <- pVals[sigs]
              consistency <- consistency[sigs]
              fdr <- fdr[sigs]

              out <- data.frame(
                  Consistency = consistency,
                  pValue = pVals,
                  FDR = fdr,
                  row.names = sigs
              )

              out <- out[order(out$pValue, out$Consistency * -1), ]

              return(out)

          })


#' Get Results of One-vs-All Differential Signature Tests
#'
#' Returns the results of running one-vs-all differential signature
#' tests for each level of every factor meta-variable.
#'
#' The 'stat' variable refers to the AUC
#'
#' The output object has a nested structure:
#'
#' List of meta-data variables -> List of variable levels -> Results Dataframe
#'
#' The results dataframe has three columns: "stat", "pValue", "FDR"
#'
#' @param object VISION object
#' @return nested list of list of data.frame (see details)
#' @export
#' @aliases getSignatureDifferential
#' @rdname getSignatureDifferential
setMethod("getSignatureDifferential", signature(object = "Vision"),
          function(object) {

              # This is almost what we want to output, but needs some massaging
              ClusterSigScores <- object@ClusterSigScores

              metaVars <- colnames(object@metaData)

              if (is.null(ClusterSigScores)){
                  return(ClusterSigScores)
              }

              to_keep <- setdiff(
                  rownames(ClusterSigScores[[1]][[1]]),
                  metaVars
                  )

              ClusterSigScores <- lapply(ClusterSigScores, function(var_res){
                  var_res <- lapply(var_res, function(var_level_res){
                      var_level_res <- var_level_res[to_keep, , drop = FALSE]
                  })
              })

              return(ClusterSigScores)
          })


#' Get Results of One-vs-All Differential Tests with Metadata Variables
#'
#' Returns the results of running one-vs-all differential
#' tests for each level of every factor meta-variable.
#'
#' For numeric meta-variables, the 'stat' is the AUC. For factor meta-variables
#' the stat is the chisq statistic comparing the two groups
#'
#' The output object has a nested structure:
#'
#' List of meta-data variables -> List of variable levels -> Results Dataframe
#'
#' The results dataframe has three columns: "stat", "pValue", "FDR"
#'
#' @param object VISION object
#' @return nested list of list of data.frame (see details)
#' @export
#' @aliases getMetaDifferential
#' @rdname getMetaDifferential
setMethod("getMetaDifferential", signature(object = "Vision"),
          function(object) {

              # This is almost what we want to output, but needs some massaging
              ClusterSigScores <- object@ClusterSigScores

              metaVars <- colnames(object@metaData)

              if (is.null(ClusterSigScores)){
                  return(ClusterSigScores)
              }

              to_keep <- intersect(
                  rownames(ClusterSigScores[[1]][[1]]),
                  metaVars
                  )

              ClusterSigScores <- lapply(ClusterSigScores, function(var_res){
                  var_res <- lapply(var_res, function(var_level_res){
                      var_level_res <- var_level_res[to_keep, , drop = FALSE]
                  })
              })

              return(ClusterSigScores)
          })
