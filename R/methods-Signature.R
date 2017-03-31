require(geometry)

setMethod("initialize", signature(.Object="Signature"),
          function(.Object, sigDict, name, source, metaData="") {
            if (missing(sigDict)) {
              stop("Missing sigDict information.")
            } else if (missing(name)) {
              stop("Missing signature name.")
            } else if (missing(source)) {
              stop("Missing source file.")
            }
            
            .Object@sigDict = sigDict
            .Object@name = name
            .Object@source = source
            .Object@metaData = metaData
            
            return(.Object);
            
          }
)

setMethod("sigEqual", signature(object="Signature"),
          function(object, compareSig) {
            if (class(compareSig) != "Signature") {
              stop("Input compareSig is not a Signature data type.")
            }
            
            if (object@name != compareSig@name) {
              return(FALSE)
            }
            
            return(object@sigDict == compareSig@sigDict)
            
          }
)

setMethod("addSigData", signature(object="Signature"),
          function(object, data) {
            object@sigDict <- list(object@sigDict, data)
            retunr(object)
          })

sigsVsProjections <- function(projections, sigScoresData, randomSigData, NEIGHBORHOOD_SIZE = 0.33) {
  #' Evaluates the significance of each signature vs. each projection. 
  #' 
  #' Paramters:
  #'  projections: ProjectionData
  #'    Maps projections to their spatial coordinates for each sample
  #'  sigScoresData: List of SignatureScores Object
  #'    Maps signature names to their value at each coordinate
  #'  randomSigData: List of SignatureScores Object
  #'    Maps randomly generated signatures to scores to be compared with the real signature scores.
  #'  
  #' Returns:
  #'  spRowLabels: List of Strings
  #'    Labels for rows of the output matrices
  #'  spColLabels: List of Strings
  #'    Labels for columns of the output matrices
  #'  sigProjMatrix: matrix (Num_Signatures x Num_Projections)
  #'    sigProj dissimilarity score
  #'  SigProjMatrixP: matrix (Num_Signatures x Num_Projections)
  
  set.seed(RANDOM_SEED)
  spRowLabels <- c()
  spRows <- c()
  spRowLabelsFactors <- c()
  spRowFactors <- c()
  spRowLabelsPNum <- c()
  spRowPNums <- c()
  
  for (sig in sigScoresData) {
    if (sig@isPrecomputed) {
      if (sig@isFactor) {
        spRowLabelsFactors <- c(spRowLabelsFactors, sig@name)
        spRowFactors <- c(spRowFactors, sig)
      } else{
        spRowLabelsPNum <- c(spRowLabelsPNum, sig@name)
        spRowPNums <- c(spRowPNums, sig)
      }
    } else {
      spRowLabels <- c(spRowLabels, sig@name)
      spRows <- c(spRows, sig)
    }
  }
  
  pKeys <- c()
  for (p in projections) {
    pKeys <- c(pKeys, p@name)
  }
  
  spColLabels <- pKeys
  spColLabels <- sort(spColLabels)
  N_SAMPLES <- length(sigScoresData[[1]]@sample_labels)
  N_SIGNATURES <- length(spRowLabels)
  N_SIGNATURE_FACTORS <- length(spRowLabelsFactors)
  N_SIGNATURE_PNUM <- length(spRowLabelsPNum)
  N_PROJECTIONS <- length(spColLabels)
  
  sigProjMatrix <- matrix(0L, nrow=N_SIGNATURES, ncol=N_PROJECTIONS)
  sigProjMatrixP <- matrix(0L, nrow=N_SIGNATURES, ncol=N_PROJECTIONS)
  
  factorSigProjMatrix <- matrix(0L, nrow=N_SIGNATURE_FACTORS, ncol=N_PROJECTIONS)
  factorSigProjMatrixP <- matrix(0L, nrow=N_SIGNATURE_FACTORS, ncol=N_PROJECTIONS)
  
  pnumSigProjMatrix <- matrix(0L, nrow=N_SIGNATURE_PNUM, ncol=N_PROJECTIONS)
  pnumSigProjMatrixP <- matrix(0L, nrow=N_SIGNATURE_PNUM, ncol=N_PROJECTIONS)
  
  # Build a matrix of all signatures
  sigScoreMatrix <- matrix(0L, nrow=N_SAMPLES, ncol=N_SIGNATURES)
  
  j <- 1
  for (sig in spRows) {
    sigScoreMatrix[,j] <- rank(sig@scores, ties.method="average")
    j <- j+1
  }

  sigScoreMatrix <- t(apply(sigScoreMatrix, 1, as.numeric))

  randomSigScoreMatrix <- matrix(0L, nrow=N_SAMPLES, ncol=length(randomSigData))
  
  j <- 1
  for (rsig in randomSigData) {
    randomSigScoreMatrix[,j] <- rank(rsig@scores, ties.method="average")
    j <- j + 1
  }
  
  randomSigScoreMatrix <- t(apply(randomSigScoreMatrix, 1, as.numeric))

  
  ### TODO: BUILD ONE HOT MATRIX FOR FACTORS 
  
  message("Evaluating signatures against projections...")
  
  i <- 1
  for (proj in projections) {
    dataLoc <- proj@pData
    distanceMatrix <- (dist(t(dataLoc), method="euclidean"))

    weights <- as.matrix(exp(-1 * ((distanceMatrix*distanceMatrix) / NEIGHBORHOOD_SIZE^2)))
    diag(weights) <- 0
    weightsNormFactor <- apply(weights, 1, sum)
    weightsNormFactor[weightsNormFactor == 0] <- 1.0
    weightsNormFactor[is.na(weightsNormFactor)] <- 1.0
    weights <-(weights / weightsNormFactor)
  
    neighborhoodPrediction <- (weights %*% sigScoreMatrix)
  
    ## Neighborhood dissimilatory score = |actual - predicted|
    dissimilarity <- abs(sigScoreMatrix - neighborhoodPrediction)
    medDissimilarity <- as.matrix(apply(dissimilarity, 2, median))
    # Calculate scores for random signatures
    randomNeighborhoodPrediction <- (weights %*% randomSigScoreMatrix)
    randomDissimilarity <- abs(randomSigScoreMatrix - randomNeighborhoodPrediction)
    randomMedDissimilarity <- as.matrix(apply(randomDissimilarity, 2, median))

    # Group by number of genes
    backgrounds <- list()
    k <- 1
    for (rsig in randomSigData) {
      numGenes <- as.character(rsig@numGenes)
      if (!(numGenes %in% names(backgrounds))) {
        backgrounds[[numGenes]] <- c(randomMedDissimilarity[k])
      } else {
        backgrounds[[numGenes]] <- c(backgrounds[[numGenes]], randomMedDissimilarity[k])
      }
      k <- k + 1
    }
   
    bgStat <- matrix(0L, nrow=length(backgrounds), ncol=3)
    k <- 1
    for (numGenes in names(backgrounds)) {
      mu_x <- mean(backgrounds[[numGenes]])
      std_x <- biasedVectorSD(as.matrix(backgrounds[[numGenes]]))
      bgStat[k, 1] <- as.numeric(numGenes)
      bgStat[k, 2] <- mu_x
      bgStat[k, 3] <- std_x
      k <- k + 1
    }
  
    
    mu <- matrix(0L, nrow=nrow(medDissimilarity), ncol=ncol(medDissimilarity))
    sigma <- matrix(0L, nrow=nrow(medDissimilarity), ncol=ncol(medDissimilarity))

    for (k in 1:length(medDissimilarity)) {
      # Find background with closest number of genes 
      numG <- spRows[[k]]@numGenes
      row_i <- which.min(abs(numG - bgStat[,1]))
      mu[k] <- bgStat[row_i, 1]
      sigma[k] <- bgStat[row_i, 2]
    }
    
    medDissimilarityPrime <- (medDissimilarity - mu) / sigma
    
    #Create CDF function for medDissmilarityPrime
    pcdf <- ecdf(medDissimilarityPrime)
    # Apply CDF function to medDissimilarityPrime pointwise
    pValues <- apply(medDissimilarityPrime, c(1, 2), pcdf)
    
    sigProjMatrix[,i] <- 1 - (medDissimilarity / N_SAMPLES)
    sigProjMatrixP[,i] <- pValues
    
   
    ## TODO 
    # Calculate significance for precomputed numerical signatures 
    # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)

    ## TODO
    # Calculate significance for Factor Signatures
    
    i <- i + 1
  }
  
  ## TODO 
  # Concatenate Factor Sig-proj entries back in 
  
  sigProjMatrixP <- apply(sigProjMatrix, 1, function(x) p.adjust(x, method="BH"))
  sigProjMatrixP[sigProjMatrixP == 0] <- 10^(-300)
  
  sigProjMatrixP <- t(apply(sigProjMatrixP, c(1,2), log10))
  colnames(sigProjMatrixP) <- spColLabels
  rownames(sigProjMatrixP) <- spRowLabels
  
  colnames(sigProjMatrix) <- spColLabels
  rownames(sigProjMatrix) <- spRowLabels
  
  return(list(spRowLabels, spColLabels, sigProjMatrix, sigProjMatrixP))
    
    
}

