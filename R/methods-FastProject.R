require(logging)
require(BiocParallel)
require(parallel)


#'  Initializes a new FastProject object.
#'
#'    @param data_file file path to expression matrix
#'    @param housekeeping file path to housekeeping data file
#'    @param signatures list of file paths to signature files (.gmt or .txt)
#'    @param scone scone normalization data file
#'    @param norm_methods normalization methods to be extracted from the scone object
#'    @param precomptued data file with precomputed signature scores (.txt)
#'    @param nofilter if TRUE, no filter applied; else filters applied. Default is FALSE
#'    @param nomodel if TRUE, no fnr curve calculated and all weights equal to 1. Else FNR and weights calculated.
#'             Default is TRUE.
#'    @param filters list of filters to compute
#'    @param debug if 1 enable debugging features, if not not. Default is 0.
#'    @param lean if TRUE run a lean simulation. Else more robust pipeline initiated. Default is FALSE
#'    @param qc if TRUE calculate QC; else not. Default is FALSE
#'    @param num_cores Number of cores to use during analysis.
#'    @param min_signature_genes Minimum number of genes required to compute a signature
#'    @param projections File containing precomputed projections for analysis
#'    @param weights Precomputed weights for each coordinate. Normally computed from the FNR curve.
#'    @param threhsold Threshold to apply for the threshold filter
#'    @param perm_wPCA If TRUE, apply permutation WPCA to calculate significant number of PCs. Else not. Default FALSE.
#'    @param sig_norm_method Method to apply to normalize the expression matrix before calculating signature scores
#'    @param sig_score_method Method to apply when calculating signature scores
#'    @param exprData expression data matrix, as opposed to specifying a data_file
#'    @param housekeepingData housekeeping data table, as opposed to specifying a housekeeping_data file
#'    @param sigData List of signatures, as opposed to specifying a list of signature files
#'    @param precomputedData List of precomputed signature scores, as opposed to specifying a precomputed file
#'    @return A FastProject object.
#'    @examples 
#'    fp <- FastProject("expression_matrix.txt", "data/Gene Name Housekeeping.txt", c("sigfile_1.gmt", "sigfile_2.txt"), precomputed="pre_sigs.txt")
setMethod("initialize", signature(.Object="FastProject"),
          function(.Object, data_file, housekeeping, signatures, scone = NULL, norm_methods = NULL, 
                   precomputed=NULL, nofilter=FALSE, nomodel=FALSE, filters=c("fano"),
                   debug=0, lean=FALSE, qc=FALSE, num_cores=1, min_signature_genes=5, projections="", 
                   weights=NULL, threshold=0, perm_wPCA=FALSE, sig_norm_method="znorm_rows", 
                   sig_score_method="weighted_avg", exprData=NULL, housekeepingData=NULL, 
                   sigData=NULL, precomputedData=NULL) {
            
            
            
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
              	.Object@data_file <- data_file
                .Object@exprData <- readExprToMatrix(data_file)
              } else {
                .Object@exprData <- exprData
              }
            }
            
            if (is.null(housekeepingData)) {
              .Object@housekeeping <- housekeeping
              .Object@housekeepingData <- readHKGToMatrix(housekeeping)
            } else {
              .Object@housekeepingData <- housekeepingData
            }
              
            if (is.null(sigData)) {
              .Object@signatures <- signatures
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
            .Object@nomodel <- nomodel
            .Object@filters <- filters
            .Object@projections <- projections
            .Object@threshold <- threshold
            .Object@sig_norm_method <- sig_norm_method
            .Object@sig_score_method <- sig_score_method
            .Object@debug = debug
            .Object@lean = lean
            .Object@perm_wPCA = perm_wPCA
            .Object@numCores = num_cores
            .Object@allData = .Object@exprData
            
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

#' Main entry point for running FastProject Analysis
#' 
#' @param object FastProject object 
#' @return FastProjectOutput object
#' @examples 
#' fp <- FastProject("expression_matrix.txt", "data/Gene Name Housekeeping.txt", c("sigfile_1.gmt", "sigfile_2.txt"), precomputed="pre_sigs.txt")
#' fpout <- Analysis(fp)
setMethod("Analyze", signature(object="FastProject"), function(object) {
  message("Beginning Analysis")
  
  ptm <- Sys.time()
  timingList <- (ptm - ptm)
  tRows <- c("Start")

  object@allData <- object@exprData

  if (ncol(object@exprData) > 35000) {	

  	  fexpr <- filterGenesFano(object@exprData)
  	  res <- applyPCA(fexpr, N=30)
  	  kn <- ball_tree_knn(t(res), 30, object@numCores)
  	  cl <- louvainCluster(kn, t(res))
  	  cl <- readjust_clusters(cl, t(res)) 

	  # Randomly select a representative cell from each cluster
  	  cell_subset <- sapply(cl, function(i) sample(i, 1))
  	  
  	  object@exprData <- object@exprData[,cell_subset]

  }

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Partition & Sample")

  # Wrap expression data frame into a ExpressionData class
  eData <- ExpressionData(object@exprData, distanceMatrix=sparseDistance)
  
  # If no filter threshold was specified, set it to 20% of samples
  if (object@threshold == 0) {
    num_samples <- ncol(getExprData(eData)) - 1
    object@threshold <- (0.2 * num_samples)
  }
  
  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Threshold")

  object@allData <- object@exprData

  filtered <- applyFilters(eData, object@threshold, object@filters)
  eData <- filtered[[1]]
  filterList <- filtered[[2]]

  originalData <- getExprData(eData)
  print(ncol(originalData))

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
  else if (all(is.na(object@weights)) || ncol(object@weights) != ncol(object@exprData)) {
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
  	if (length(s) != ncol(object@exprData)) {
		s@scores <- s@scores[colnames(object@exprData)]
		s@sample_labels <- colnames(object@exprData)
	}
    sigScores <- c(sigScores, s)
    sigList<- c(sigList, Signature(list(), s@name, "", "", isPrecomputed=TRUE, isFactor=s@isFactor, cluster=0))
    sigNames <- c(sigNames, s@name)
  }
  names(sigList) <- sigNames
  names(sigScores) <- sigNames

  # Remove any signatures that didn't compute correctly
  toRemove <- lapply(sigScores, function(x) all(x@scores == c(-Inf)))
  sigScores <- sigScores[names(which(toRemove==F))]

  ## Convert Sig Scores to matrix
  names <- c()
  sigMatrix <- matrix(0L, nrow=length(sigScores), ncol=length(sigScores[[1]]@scores))
  for (sig in 1:length(sigScores)) {
    names <- c(names, sigScores[[sig]]@name)
    scores <- t(as.matrix(sigScores[[sig]]@scores))
    sigMatrix[sig,] <- scores
  }

  rownames(sigMatrix) <- names
  colnames(sigMatrix) <- colnames(object@exprData)
  
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

    projectData <- generateProjections(eData, object@weights, filter, inputProjections <- c(), lean=object@lean, perm_wPCA = object@perm_wPCA, numCores = object@numCores)
    projs <- projectData[[1]]
    g <- projectData[[2]]
    PPT <- projectData[[3]]
    pca_res <- projectData[[4]]
    loadings <- projectData[[5]]
    
    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("Pr. ", filter))

    message("Computing significance of signatures...")
    sigVProj <- sigsVsProjections(projs, sigScores, randomSigScores)
    
    sigKeys <- sigVProj[[1]]
    projKeys <- sigVProj[[2]]
    sigProjMatrix <- sigVProj[[3]]
    pVals <- sigVProj[[4]]
    
    timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
    tRows <- c(tRows, paste0("SigVPr ", filter))

	pearsonCorr <- lapply(1:nrow(sigMatrix), function(i) {
		lapply(1:nrow(pca_res), function(j) {
			return(calcPearsonCorrelation(sigMatrix[i,], pca_res[j,]))
		})
	})

	pearsonCorr <- lapply(pearsonCorr, unlist)
	pearsonCorr <- matrix(unlist(pearsonCorr), nrow=nrow(sigMatrix), ncol=nrow(pca_res), byrow=T)

	rownames(pearsonCorr) <- rownames(sigMatrix)
	colnames(pearsonCorr) <- rownames(pca_res)

	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
	tRows <- c(tRows, paste("Pearson Correlation", filter))

    projData <- ProjectionData(filter=filter, projections=projs, genes=g, keys=projKeys, 
                                          sigProjMatrix=sigProjMatrix, pMatrix=pVals, PPT=PPT, fullPCA=pca_res,
                                          pearsonCorr=pearsonCorr, loadings=loadings)

    projDataList[[filter]] <- projData
      
  }
    

  sigList <- sigList[rownames(sigMatrix)]
  
  message("Clustering Signatures...")
  sigClusters <- clusterSignatures(sigList, sigMatrix, k=10)
    

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "ClusterSignatures")


  slots <- slotNames("FastProject")
  fpParams <- list()
  for (s in slots) {
	fpParams[[s]] <- slot(object, s)
  }

  fpOut <- FastProjectOutput(eData, projDataList, sigMatrix, sigList, sigClusters, fpParams)

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Final")

  timingList <- as.matrix(timingList)
  rownames(timingList) <- tRows

  #return(list(fpOut, timingList))
  return(fpOut)
}
)

#' Creates a new FastProject object with the new expression matrix. Called from Server.R.
#' 
#' @param fpParams List of FastProject parameters
#' @param nexpr New expression matrix for the new FastProject object
#' @return a new FastProject object
#' @examples 
#'   fp <- FastProject("expression_matrix.txt", "data/Gene Name Housekeeping.txt", c("sigfile_1.gmt", "sigfile_2.txt"), precomputed="pre_sigs.txt")
#'   slots <- slotNames(object)
#'   fpParams <- list()
#'   for (s in slots) {
#'     fpParams <- c(fpParams, slot(object, s))
#'   }
#'   subset_cells <- sample(1:ncol(fp@exprData), 200)
#'   nexpr <- fp@exprData[,subset_cells]
#'   names(fpParams) <- slots
#'   nfp <- extraAnalysisFastProject(fpParams, nexpr)

createNewFP <- function(fpParams, nexpr) {
	nfp <- FastProject(exprData=nexpr, sigData=fpParams[["sigData"]],
						housekeepingData=fpParams[["housekeepingData"]])

	slots <- names(fpParams)	
	for (s in slots) {
		if (s != "exprData" && s != "sigData" && s != "housekeepingData") {
			slot(nfp, s) <- fpParams[[s]]
		}
	}

	return(nfp)
}
  
