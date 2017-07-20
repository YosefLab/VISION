require(logging)
require(BiocParallel)
require(parallel)

setMethod("initialize", signature(.Object="FastProject"),
          function(.Object, data_file, housekeeping, signatures, scone = NULL, norm_methods = NULL, 
                   precomputed=NULL, output_dir = "FastProject_Output", nofilter=FALSE, nomodel=FALSE, filters=c("fano"),
                   all_sigs=FALSE, debug=0, lean=TRUE, subsample_size=0, qc=FALSE, num_cores=1, approximate=F, optClust=0, 
                   min_signature_genes=5, projections="", weights=NULL, threshold=0, perm_wPCA=FALSE,
                   sig_norm_method="znorm_rows", sig_score_method="weighted_avg", exprData=NULL, 
                   housekeepingData=NULL, sigData=NULL, precomputedData=NULL) {
            
            ## Make sure that the minimum files are being read in.
            if (missing(data_file) && is.null(exprData) && is.null(scone)) {
              stop("Missing expression data file or SCONE object in input.")
            } else if (missing(housekeeping) && is.null(housekeepingData)) {
              stop("Missing housekeeping data file in input.")
            } else if (missing(signatures) && is.null(sigData)) {
              stop("Missing signature data file in input.") 
            }
            
            if (is.null(scone)) {
              if (is.null(exprData)) {
                .Object@exprData <- readExprToMatrix(data_file)
              } else {
                .Object@exprData <- exprData
              }
            }
            
            if (is.null(housekeepingData)) {
              .Object@housekeepingData <- readHKGToMatrix(housekeeping)
            } else {
              .Object@housekeepingData <- housekeepingData
            }
              
            if (is.null(sigData)) {
              .Object@sigData <- readSignaturesInput(signatures)
            } else {
              .Object@sigData <- sigData
            }
            
            if (is.null(weights)) {
              .Object@weights <- matrix(NA, nrow=10, ncol=0)
            } else {
              .Object@weights <- weights
            }
            
            if (!is.null(precomputed)) {
              .Object@precomputedData <- readPrecomputed(precomputed, colnames(.Object@exprData))
            }
            
            .Object@nofilter <- nofilter
            .Object@output_dir <- output_dir
            .Object@nomodel <- nomodel
            .Object@subsample_size <- subsample_size
            .Object@filters <- filters
            .Object@projections <- projections
            .Object@threshold <- threshold
            .Object@sig_norm_method <- sig_norm_method
            .Object@sig_score_method <- sig_score_method
            .Object@debug = debug
            .Object@lean = lean
            .Object@perm_wPCA = perm_wPCA
            .Object@approximate = approximate
            .Object@numCores = num_cores
            .Object@optClust = optClust
            
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
  
  ptm <- Sys.time()
  timingList <- (ptm - ptm)
  tRows <- c("Start")

  # Wrap expression data frame into a ExpressionData class
  eData <- ExpressionData(object@exprData, distanceMatrix=sparseDistance)
  
  holdouts <- NULL
  if (object@subsample_size != 0) {
    split_samples <- splitSamples(getExprData(eData), object@subsample_size)
    holdouts <- split_samples[[1]]
    subdata <- split_samples[[2]]
    eData <- updateExprData(eData, split_samples[[2]])
  }
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Holdouts")

  # If no filter threshold was specified, set it to 20% of samples
  if (object@threshold == 0) {
    num_samples <- ncol(getExprData(eData)) - 1
    object@threshold <- (0.2 * num_samples)
  }
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Threshold")

  originalData <- getExprData(eData)

  filtered <- applyFilters(eData, object@threshold, object@filters)
  eData <- filtered[[1]]
  filterList <- filtered[[2]]
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Filter") 

  falseneg_out <- createFalseNegativeMap(originalData, object@housekeepingData, object@debug)
  func <- falseneg_out[[1]]
  params <- falseneg_out[[2]]
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "FN Function")

  if (object@nomodel) {
    object@weights <- matrix(1L, nrow=nrow(originalData), ncol=ncol(originalData))
    rownames(object@weights) <- rownames(originalData)
    colnames(object@weights) <- colnames(originalData)
  }
  else if (all(is.na(object@weights))) {
    message("Computing weights from False Negative Function...")
    object@weights <- computeWeights(func, params, eData)
  }

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Weights")

  zero_locations <- which(getExprData(eData) == 0.0, arr.ind=TRUE) 
    
  normalizedData <- getNormalizedCopy(eData, object@sig_norm_method)
  eData <- updateExprData(eData, normalizedData)
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Normalize")    

  ## Define single signature evaluation for lapply method
  singleSigEval <- function(s) {
	x <- c(-Inf)
	if (object@sig_score_method=="naive") {
		tryCatch({
			x <- naiveEvalSignature(eData, s, object@weights, object@min_signature_genes)
		}, error=function(e){})
	} else if (object@sig_score_method=="weighted_avg") {
		tryCatch({
			x <- weightedEvalSignature(eData, s, object@weights, object@min_signature_genes)
		}, error=function(e){})
	}
	return(x)
  }	

  # Score user defined signatures with defined method (default = weighted)
  sigScores <- c()
  if (object@sig_score_method == "naive") {
    message("Applying naive signature scoring method...")
  } else if (object@sig_score_method == "weighted_avg") {
    message("Applying weighted signature scoring method...")
  }

  sigScores <- bplapply(object@sigData, singleSigEval, BPPARAM=MulticoreParam(workers=object@numCores))
  
  sigList <- object@sigData
  sigNames <- names(object@sigData)
  for (s in object@precomputedData) {
    sigScores <- c(sigScores, s)
    sigList<- c(sigList, Signature(list(), s@name, "", "", isPrecomputed=TRUE, isFactor=s@isFactor, cluster=0))
    sigNames <- c(sigNames, s@name)
  }
  names(sigList) <- sigNames
  names(sigScores) <- sigNames

  # Remove any signatures that didn't compute correctly
  toRemove <- lapply(sigScores, function(x) all(x@scores == c(-Inf)))
  sigScores <- sigScores[names(which(toRemove==F))]
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Sig Scores")

  # Construct random signatures for background distribution
  randomSigs <- c()
  randomSizes <- c(5, 10, 20, 50, 100, 200)
  for (size in randomSizes) {
    message("Creating random signature of size ", size, "...")
    for (j in 1:2000) {
      newSigGenes <- sample(rownames(getExprData(eData)), size)
      newSigSigns <- rep(1, size)
      names(newSigSigns) <- newSigGenes
      newSig <- Signature(newSigSigns, paste0("RANDOM_BG_", size, "_", j), 'x')
      randomSigs <- c(randomSigs, newSig)
    }
  }
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Random Sigs")

  # Compute signature scores for random signatures generated
  randomSigScores <- c()
  if (object@sig_score_method == "naive") {
    message("Applying naive signature scoring method...")
  } else if (object@sig_score_method == "weighted_avg") {
    message("Applying weighted signature scoring method...")
  }
	
  randomSigScores <- bplapply(randomSigs, singleSigEval,BPPARAM=MulticoreParam(workers=object@numCores))
  names(randomSigScores) <- names(randomSigs)
  

  # Remove random signatures that didn't compute correctly
  toRemove <- lapply(randomSigScores, function(x) all(x@scores == c(-Inf)))
  randomSigScores <- randomSigScores[which(toRemove==F)]
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Rand Sig Scores")

  # Apply projections to filtered gene sets, create new projectionData object
  projDataList <- list()
  for (filter in filterList) {
    message("Filter level: ", filter)

    message("Projecting data into 2 dimensions...")

    projectData <- generateProjections(eData, object@weights, filter, inputProjections <- c(), lean=object@lean, perm_wPCA = object@perm_wPCA, numCores = object@numCores, optClust = object@optClust, approximate = object@approximate)
    projs <- projectData[[1]]
    g <- projectData[[2]]
    PPT <- projectData[[3]]
    
    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("Pr. ", filter))

    message("Computing significance of signatures...")
    #return(list(projs, sigScores, randomSigScores))
    sigVProj <- sigsVsProjections(projs, sigScores, randomSigScores)
    
    sigKeys <- sigVProj[[1]]
    projKeys <- sigVProj[[2]]
    sigProjMatrix <- sigVProj[[3]]
    pVals <- sigVProj[[4]]
    
    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("SigVPr ", filter))

    projData <- ProjectionData(filter=filter, projections=projs, genes=g, keys=projKeys, 
                                          sigProjMatrix=sigProjMatrix, pMatrix=pVals, PPT=PPT)

    projDataList[[filter]] <- projData
      
  }
    
  ## Convert Sig Scores to matrix
  names <- c()
  sigMatrix <- matrix(0L, nrow=length(sigScores), ncol=length(sigScores[[1]]@scores))
  for (sig in 1:length(sigScores)) {
    names <- c(names, sigScores[[sig]]@name)
    scores <- t(as.matrix(sigScores[[sig]]@scores))
    sigMatrix[sig,] <- scores
  }
  
  rownames(sigMatrix) <- names
  colnames(sigMatrix) <- colnames(fp@exprData)

  sigList <- sigList[rownames(sigMatrix)]
  
  sigClusters <- clusterSignatures(sigList, sigMatrix, k=10)
    

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "ClusterSignatures")

  fpOut <- FastProjectOutput(eData, projDataList, sigMatrix, sigList, sigClusters)

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Final")

  timingList <- as.matrix(timingList)
  rownames(timingList) <- tRows

  #return(list(fpOut, timingList))
  return(fpOut)
}
  
)
