require(logging)

setMethod("initialize", signature(.Object="FastProject"),
          function(.Object, data_file, housekeeping, signatures, precomputed=NULL,
                   output_dir = "FastProject_Output", nofilter=FALSE, nomodel=FALSE, pca_filter=FALSE,
                   all_sigs=FALSE, debug=0, lean=FALSE, subsample_size=0, qc=FALSE, 
                   min_signature_genes=5, projections="", weights=NULL, threshold=0, 
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
              .Object@exprData <- readExprToMatrix(data_file)
            } else {
              .Object@exprData <- exprData
            }
            
            if (is.null(housekeepingData)) {
              .Object@housekeepingData <- readHKGToMatrix(housekeeping)
            } else {
              .Object@housekeepingData <- housekeepingData
            }
              
            if (is.null(sigData)) {
              .Object@sigData <- readSignaturesGmtToMatrix(signatures)
            } else {
              .Object@sigData <- sigData
            }
            
            if (is.null(weights)) {
              .Object@weights <- matrix(NA, nrow=10, ncol=0)
            } else {
              .Object@weights <- weights
            }
            
            .Object@nofilter <- nofilter
            .Object@output_dir <- output_dir
            .Object@nomodel <- nomodel
            .Object@subsample_size <- subsample_size
            .Object@pca_filter <- pca_filter
            .Object@projections <- projections
            .Object@threshold <- threshold
            .Object@sig_norm_method <- sig_norm_method
            .Object@sig_score_method <- sig_score_method
            .Object@debug = debug
            
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

  holdouts <- NULL
  if (object@subsample_size != 0) {
    split_samples <- splitSamples(getExprData(eData), object@subsample_size)
    holdouts <- split_samples[[1]]
    subdata <- split_samples[[2]]
    eData <- updateExprData(eData, split_samples[[2]])
  }
  
  # If no filter threshold was specified, set it to 20% of samples
  if (object@threshold == 0) {
    num_samples <- ncol(getExprData(eData)) - 1
    object@threshold <- (0.2 * num_samples)
  }
  
  originalData <- copy(getExprData(eData))

  filtered <- applyFilters(eData, object@threshold, object@nofilter, object@lean)
  eData <- filtered[[1]]
  filterList <- filtered[[2]]
  
  if (!object@nomodel) {
    
    falseneg_out <- createFalseNegativeMap(originalData, object@housekeepingData, object@debug)
    func <- falseneg_out[[1]]
    params <- falseneg_out[[2]]
    
    if (all(is.na(object@weights))) {
      object@weights <- computeWeights(func, params, eData)
    }
    
    zero_locations <- which(getExprData(eData) == 0.0, arr.ind=TRUE) 
    
    normalizedData <- getNormalizedCopy(eData, object@sig_norm_method)
    eData <- updateExprData(eData, normalizedData)
    
    
    # Score user defined signatures with defined method (default = weighted)
    sig_scores <- c()
    if (object@sig_score_method == "naive") {
      message("Applying naive signature scoring method...")
      for (sig in fp@sigData) {
        tryCatch({
          sig_scores <- c(sig_scores, naiveEvalSignature(eData, 
                                          sig, fp@weights, object@min_signature_genes))
        }, error=function(e){})
      }
    } else if (object@sig_score_method == "weighted_avg") {
      message("Applying weighted signature scoring method...")
      for (sig in fp@sigData) {
        tryCatch({
          sig_scores <- c(sig_scores, weightedEvalSignature(eData, 
                                          sig, fp@weights, object@min_signature_genes))
        }, error=function(e){})
      }
    }
    
    # Construct random signatures for background distribution
    randomSigs <- c()
    randomSizes <- c(5, 10, 20, 50, 100, 200)
    for (size in randomSizes) {
      message("Creating random signature of size ", size, "...")
      for (j in 1:3000) {
        newSigGenes <- sample(rownames(getExprData(eData)), size)
        newSigExpr <- rep(1, size)
        names(newSigExpr) <- newSigGenes
        newSig <- Signature(newSigExpr, paste0("RANDOM_BG_", size, "_", j), 'x')
        randomSigs <- c(randomSigs, newSig)
      }
    }
    
    
    # Compute signature scores for random signatures generated
    randomSigScores <- c()
    if (object@sig_score_method == "naive") {
      message("Applying naive signature scoring method on random signatures...")
      for (sig in randomSigs) {
        tryCatch({
          randomSigScores <- c(randomSigScores, naiveEvalSignature(eData, 
                                                         sig, fp@weights, object@min_signature_genes))
        }, error=function(e){})
      }
    } else if (object@sig_score_method == "weighted_avg") {
      message("Applying weighted signature scoring method on random signatures...")
      for (sig in randomSigs) {
        tryCatch({
          randomSigScores <- c(randomSigScores, weightedEvalSignature(eData, 
                                                            sig, fp@weights, object@min_signature_genes))
        }, error=function(e){})
      }
    }
  
    
    # Apply projections to filtered gene sets, create new projectionData object
    for (filter in filterList) {
      message("Filter level: ", filter)
      
      message("Projecting data into 2 dimensions...")
      
      projectData <- generateProjections(edata, filter)
      projections <- projectData[[1]]
      genes <- projectData[[2]]
      
      clusters <- defineClusters(projections)
      
      pKeys <- c()
      for (p in projections) {
        pKeys <- c(pKeys, p@name)
      }
      
      projData <- projectionData(filter=filter, pca=TRUE, projections=projections, genes=genes, keys=pKeys)
      
      
    }
    
    
    return(getExprData(eData))
  }
  
  }
)




