require(logging)

setMethod("initialize", signature(.Object="FastProject"),
          function(.Object, data_file, housekeeping, signatures, precomputed,
                   output_dir = "FastProject_Output", nofilter = FALSE, nomodel = FALSE, pca_filter = FALSE,
                   all_sigs = FALSE, debug = FALSE, lean = FALSE, subsample_size = 1000, min_signature_genes = 5,
                   projections = NULL, weights = NULL, threshold = 0, sig_norm_method = "znorm_rows",
                   sig_score_method ="weighted avg", exprData = NULL, housekeepingData = NULL,
                   sigData = NULL) {
            
            ## Make sure that the minimum files are being read in.
            if (missing(data_file) && is.null(exprData)) {
              stop("Missing expression data file in input.")
            } else if (missing(housekeeping) && is.null(housekeepingData)) {
              stop("Missing housekeeping data file in input.")
            } else if (missing(signatures) && is.null(sigData)) {
              stop("Missing signature data file in input.")
            }
            
            if (is.null(exprData)) {
              exprData <- readTextToMatrix(data_file)
            }
            if (is.null(housekeepingData)) {
              housekeepingData = readTextToMatrix(housekeeping)
            }
            if (is.null(sigData)) {
              sigData = readSignaturesGmtToMatrix(signatures)
            }
            
            createOutputDirectory(.Object)
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
  return(NULL)
  
  }
)




