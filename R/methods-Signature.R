require(geometry)

setMethod("initialize", signature(.Object="Signature"),
          function(.Object, sigDict, name, source, metaData="", isPrecomputed=F, isFactor=F) {
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
            .Object@isPrecomputed = isPrecomputed
            .Object@isFactor = isFactor
            
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


BG_DIST <- matrix(0L, nrow=0, ncol=0)

getBGDist <- function(N_SAMPLES, NUM_REPLICATES) {
  
  if (nrow(BG_DIST) != N_SAMPLES || ncol(BG_DIST) != NUM_REPLICATES) {
    set.seed(RANDOM_SEED)
    BG_DIST <- matrix(rnorm(N_SAMPLES*NUM_REPLICATES), nrow=N_SAMPLES, ncol=NUM_REPLICATES)
    BG_DIST <- apply(BG_DIST, 2, order)
  }
  
  return(BG_DIST)
  
}

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
  #'  SigProjMatrix_P: matrix (Num_Signatures x Num_Projections)
  
  set.seed(RANDOM_SEED)
  spRowLabels <- c()
  spRows <- c()
  spRowLabelsFactors <- c()
  spRowFactors <- c()
  spRowLabelsPNum <- c()
  spRowPNums <- c()
  
  precomputedSigs <- c()
  precomputedNumerical <- c()
  precomputedFactor <- c()
  for (sig in sigScoresData) {
    if (sig@isPrecomputed) {
      if (sig@isFactor) {
        spRowLabelsFactors <- c(spRowLabelsFactors, sig@name)
        spRowFactors <- c(spRowFactors, sig)
        precomputedSigs <- c(precomputedSigs, sig)
        precomputedFactor <- c(precomputedFactor, sig)
      } else{
        spRowLabelsPNum <- c(spRowLabelsPNum, sig@name)
        spRowPNums <- c(spRowPNums, sig)
        precomputedSigs <- c(precomputedSigs, sig)
        precomputedNumerical <- c(precomputedNumerical, sig)
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
  sigProjMatrix_P <- matrix(0L, nrow=N_SIGNATURES, ncol=N_PROJECTIONS)
  
  factorSigProjMatrix <- matrix(0L, nrow=N_SIGNATURE_FACTORS, ncol=N_PROJECTIONS)
  factorSigProjMatrix_P <- matrix(0L, nrow=N_SIGNATURE_FACTORS, ncol=N_PROJECTIONS)
  
  pnumSigProjMatrix <- matrix(0L, nrow=N_SIGNATURE_PNUM, ncol=N_PROJECTIONS)
  pnumSigProjMatrix_P <- matrix(0L, nrow=N_SIGNATURE_PNUM, ncol=N_PROJECTIONS)
  
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

  
<<<<<<< HEAD
  ### BUILD ONE HOT MATRIX FOR FACTORS 
=======
  ### TODO: BUILD ONE HOT MATRIX FOR FACTORS 
>>>>>>> 57e2cb75e66c2f61c1d8afea5d8fcb1c211c8e46
  factorSigs <- list()
  for (s in precomputedFactor) {
    fValues <- s@scores
    fLevels <- unique(fValues)
    factorFreq <- matrix(0L, ncol=length(fLevels))
    factorMatrix <- matrix(0L, nrow=N_SAMPLES, ncol=length(fLevels))
    j <- 1
    for (fval in fLevels) {
      factorMatrixRow <- matrix(0L, nrow=N_SAMPLES, ncol=1)
      equal_ii <- which(fValues == fval)
      factorMatrixRow[equal_ii] <- 1
      factorFreq[j] <- length(equal_ii) / length(fValues)
      factorMatrix[,j] <- factorMatrixRow
      j <- j+1 
    }

    factorSigs[[s@name]] <- list(fLevels, factorFreq, factorMatrix)
  
  }
  
  
  
  message("Evaluating signatures against projections...")
  
  i <- 1
  for (proj in projections) {
    dataLoc <- proj@pData
    distanceMatrix <- (dist(t(dataLoc), method="euclidean"))

    weights <- as.matrix(exp(-1 * (distanceMatrix*distanceMatrix) / NEIGHBORHOOD_SIZE^2))
    diag(weights) <- 0
    weightsNormFactor <- apply(weights, 1, sum)
    weightsNormFactor[weightsNormFactor == 0] <- 1.0
    weightsNormFactor[is.na(weightsNormFactor)] <- 1.0
    weights <- weights / weightsNormFactor
    print(proj@name)
    print(dim(sigScoreMatrix))
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
      mu[k] <- bgStat[row_i, 2]
      sigma[k] <- bgStat[row_i, 3]
    }

    #Create CDF function for medDissmilarityPrime and apply CDF function to medDissimilarityPrime pointwise
    pValues <- pnorm( ((medDissimilarity - mu) / sigma))

    sigProjMatrix[,i] <- 1 - (medDissimilarity / N_SAMPLES)
    sigProjMatrix_P[,i] <- pValues
   
    # Calculate significance for precomputed numerical signatures 
    # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)
    j <- 1
    for (s in precomputedNumerical) {
      
      sigScores <- rank(s@scores)
      if (all(sigScores == sigScores[1])) {
        pnumSigProjMatrix[j,i] <- 0.0
        pnumSigProjMatrix_P[j,i] <- 1.0
        next
      }
      
      sigPredictions <- (weights %*% sigScores)
      dissimilarity <- abs(sigScores - sigPredictions)
      medDissimilarity <- as.matrix(apply(dissimilarity, 2, median))
      
      #Compute a background for numerical signatures
      NUM_REPLICATES <- 10000
      randomSigValues <- get_bg_dist(N_SAMPLES, NUM_REPLICATES)
      bgValues <- as.vector(t(sigScores))[randomSigValues]
      randomPredictions <- (weights %*% bgValues)
      rDissimilarity <- abs(bgValues <- randomPredictions)
      randomScores <- as.matrix(apply(rDissimilarity, 2, median))
      
      mu <- mean(randomScores)
      sigma <- biasedVectorSD(random_scores)
      if (sigma != 0) {
        p_value <- pnorm((medDissimilarity - mu) / sigma)
      } else {
        p_value <- 1.0
      }
      
      pnumSigProjMatrix[j,i] <- 1 - (medDissimilarity / N_SAMPLES)
      pnumSigProjMatrix_P[j,i] <- p_value
      
      j <- j + 1
    }
    
    
    # Calculate signficance for Factor signatures
    j <- 1
    for (s in factorSigs) {
      fLevels <- s[[1]]
      factorFreq <- s[[2]]
      fMatrix <- s[[3]]
      
      if (1 %in% factorFreq) {
        factorSigProjMatrix[j,i] <- 0.0
        factorSigProjMatrix_P[j,i] <- 1.0
        next
      }
      
      N_LEVELS <- length(fLevels)
      factorPredictions <- (weights %*% fMatrix)
      
      dissimilarity <- 1 - as.matrix(apply(fMatrix %*% t(factorPredictions), 1, sum))

      medDissimilarity <- as.matrix(apply(dissimilarity, 2, median))
      
      NUM_REPLICATES <- 1000
      columnSamples <- sample(N_LEVELS, NUM_REPLICATES, replace=TRUE)
      columnAssignments <- factorFreq[columnSamples]
      dim(columnAssignments) <- c(1, NUM_REPLICATES)
      randFactors <- matrix(runif(N_SAMPLES*NUM_REPLICATES), nrow=N_SAMPLES, ncol=NUM_REPLICATES)
      randPredictions <- (weights %*% randFactors)
      for (ii in 1:nrow(randPredictions)) {
        randPredictions[ii,] <- sample(randPredictions[ii,])
      }
      randMedDissimilarity <- as.matrix(apply(1-randPredictions, 2, median))
      
      mu <- mean(randMedDissimilarity)
      sigma <- biasedVectorSD(randMedDissimilarity)
      if (sigma == 0) {
        p_value <- 1.0
      } else{
        p_value <- pnorm((medDissimilarity - mu) / sigma)
      }
      
      factorSigProjMatrix[j,i] <- 1 - medDissimilarity
      factorSigProjMatrix_P[j,i] <- p_value
    }
    
    i <- i + 1
  }
  

  # Concatenate Factor Sig-proj entries back in 
  sigProjMatrix <- rbind(sigProjMatrix, factorSigProjMatrix, pnumSigProjMatrix)
  sigProjMatrix_P <- rbind(sigProjMatrix_P, factorSigProjMatrix_P, pnumSigProjMatrix_P)
  spRowLabels <- c(spRowLabels, spRowLabelsFactors, spRowLabelsPNum)
  
  original_shape <- dim(sigProjMatrix_P)
  sigProjMatrix_P <- p.adjust(sigProjMatrix_P, method="BH")
  dim(sigProjMatrix_P) <- original_shape
  sigProjMatrix_P[sigProjMatrix_P == 0] <- 10^(-300)
  
  sigProjMatrix_P <- as.matrix(log10(sigProjMatrix_P))
  
  colnames(sigProjMatrix_P) <- spColLabels
  rownames(sigProjMatrix_P) <- spRowLabels
  
  colnames(sigProjMatrix) <- spColLabels
  rownames(sigProjMatrix) <- spRowLabels
  
  return(list(spRowLabels, spColLabels, sigProjMatrix, sigProjMatrix_P))
}

