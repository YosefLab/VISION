require(logging)

setMethod("initialize", signature(.Object="FastProject"),
          function(.Object, data_file, housekeeping, signatures, precomputed=NULL,
                   output_dir = "FastProject_Output", nofilter=FALSE, nomodel=FALSE, pca_filter=FALSE,
                   all_sigs=FALSE, debug=FALSE, lean=FALSE, subsample_size=FALSE, qc=FALSE, 
                   min_signature_genes=5, projections="", weights="", threshold=0, 
                   sig_norm_method="znorm_rows", sig_score_method="weighted_avg", exprData=NULL, 
                   housekeepingData=NULL, sigData=NULL) {
            
            ## Make sure that the minimum files are being read in.
            if (missing(data_file) && is.null(exprData)) {
              stop("Missing expression data file in input.")
            } else if (missing(housekeeping) && is.null(housekeepingData)) {
              stop("Missing housekeeping data file in input.")
            } else if (missing(signatures) && is.null(sigData)) {
              stop("Missing signature data file in input.")
            }
            
            if (is.null(exprData)) {
              .Object@exprData <- readTextToMatrix(data_file)
            }
            if (is.null(housekeepingData)) {
              .Object@housekeepingData = readTextToMatrix(housekeeping)
            }
            if (is.null(sigData)) {
                .Object@sigData = readSignaturesGmtToMatrix(signatures)
              }
          
            .Object@nofilter <- nofilter
            .Object@output_dir <- output_dir
            .Object@nomodel <- nomodel
            .Object@pca_filter <- pca_filter
            .Object@projections <- projections
            .Object@weights <- weights
            .Object@threshold <- threshold
            .Object@sig_norm_method <- sig_norm_method
            .Object@sig_score_method <- sig_score_method
            
            #createOutputDirectory(.Object)
            return(.Object) 
          }
)

setMethod("createOutputDirectory", "FastProject", function(object) {
  #' Creates the output directory structure
  mainDir <- getwd()
  if (dir.exists(file.path(mainDir, object@output_dir))) {
    i = 1
    while(TRUE) {
      dir_name <- paste0(object@output_dir, toString(i))
      if (!dir.exists(file.path(mainDir, dir_name))) {
        break
      } else {
        i = i + 1
      }    
    }
  } else {
    dir_name = object@output_dir
  }
  
  dir.create(file.path(mainDir, dir_name))
  
  logger <- getLogger("FastProject")
  
  #logger$setLevel(logging.INFO)
  #fh <- FileHandler(file.path(mainDir, dir_name, 'fastproject.log'))
  #fh$setFormatter(Formatter('%(asctime)s %(message)s'))
  #logger.addHandler(fh)
  
  #for (x in slotNames(object)) {
  #  logger$info(paste0(x, toString(object@x)))
  #}
}
)

setMethod("Analyze", signature(object="FastProject"), function(object) {
  
  # Wrap expression data frame into a ExpressionData class
  eData <- ExpressionData(object@exprData)
  
  # TODO: ADD SUBSAMPLE CLASS, THEN ADD THRESHOLD CALCULATION
  #holdouts <- NULL
  #split_samples <- SubSample.split_samples(exprData, .Object@subsample_size)
  
  # If no filter threshold was specified, set it to 20% of samples
  #if (.Object@threshold == 0) {
  #  
  #}
  
  originalData <- copy(getExprData(eData))
  
  applyFilters(eData, object@threshold, object@nofilter, object@lean)
  
  if (!object@nomodel) {
    
    falseneg_out <- createFalseNegativeMap(originalData, object@housekeepingData)
    
    func <- falseneg_out[1]
    params <- falseneg_out[2]
    
    norm_methods = c("none" = noNormalization, 
                     "znorm_columns" = colNormalization,
                     "znorm_rows" = rowNormalization,
                     "znorm_rows_then_columns" = rowAndColNormalization,
                     "rank_norm_columns" = colRankNormalization)
    
    normalizedData <- getNormalizedCopy(norm_methods[.Object@sig_norm_method])
    
  }
  
  
  
  
  }
)




