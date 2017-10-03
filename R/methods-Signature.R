require(geometry)
require(pROC)
require(entropy)
require(Hmisc)
require(mclust)

#' Initialize a new Signature object.
#'
#' @param sigDict Named list of signs for each gene in the signature
#' @param name Name of the signature
#' @param source File from which this signature was read from
#' @param metaData Metadata pertinent to signature
#' @param isPreocomputed If TRUE indicates that this signature was precomputed. Else not precomputed. Default is F.
#' @param isFactor If TRUE indicates that htis signature is a Factor, else not a factor. Default is F.
#' @param cluster Number representing which cluster this signature is a part of. Default is 0.
#' @return Signature object
setMethod("initialize", signature(.Object="Signature"),
          function(.Object, sigDict, name, source, metaData="", isPrecomputed=F, isFactor=F, cluster=0) {
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
            .Object@cluster = cluster


            return(.Object);

          }
)

#' Method for testing whether or not signatures are equal to one another (if all signs, genes and name are the same.)
#'
#' @param object Signature object
#' @param compareSig Other signature to compare against
#' @return TRUE if sigs are equal, else FALSE
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

#' Add signature data to the sigDict of this Signature Object
#'
#' @param object Signature object
#' @param data Additional signature data, in the form of a named list
#' @return Updated Signature object.
setMethod("addSigData", signature(object="Signature"),
          function(object, data) {
            object@sigDict <- list(object@sigDict, data)
            return(object)
          })


BG_DIST <- matrix(0L, nrow=0, ncol=0)

#' Generates, or if already generated retrieves, a random Background Distribution
#'
#' @param N_SAMPLES Number of samples to generate the distribution for
#' @param NUM_REPLICATES Number of replicates to generate for the background distribution
#' @return Random matrix with dimensions N_SAMPLES x NUM_REPLICATES with row ordered in ascending order
#' @examples
#' bgdist <- getBGDist(430, 1000)
getBGDist <- function(N_SAMPLES, NUM_REPLICATES) {

  if (nrow(BG_DIST) != N_SAMPLES || ncol(BG_DIST) != NUM_REPLICATES) {
    set.seed(RANDOM_SEED)
    BG_DIST <- matrix(rnorm(N_SAMPLES*NUM_REPLICATES), nrow=N_SAMPLES, ncol=NUM_REPLICATES)
    BG_DIST <- apply(BG_DIST, 2, order)
  }

  return(BG_DIST)

}

#' Evaluates the significance of each signature vs. each projection.
#'
#'  @param projections Maps projections to their spatial coordinates for each sample
#'  @param sigScoresData: List of SignatureScores Object, mapping signature names to their value at each coordinate
#'  @param randomSigData: List of SignatureScores Object, mapping randomly generated signatures to scores to be compared with the real signature scores.
#'  @return Labels for rows of the output matrices
#'  @return Labels for columns of the output matrices
#'  @return Matrix (Num_Signatures x Num_Projections) sigProj dissimilarity score
#'  @return SigProjMatrix_P matrix (Num_Signatures x Num_Projections) sigProj p values
sigsVsProjections <- function(projections, sigScoresData, randomSigData, NEIGHBORHOOD_SIZE = 0.33, BPPARAM=bpparam()) {
  ptm <- Sys.time()
  tRows <- c()

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

  timingList <- c(difftime(Sys.time(), ptm, units="secs"))
  tRows <- c(tRows, "Initialize all Matrices")

  # Build a matrix of all signatures
  sigScoreMatrix <- matrix(unlist(bplapply(spRows, function(sig) { rank(sig@scores, ties.method="average") },
  										   BPPARAM=BPPARAM)), nrow=N_SAMPLES, ncol=length(spRows))


  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "SigScoreMatrix")


  randomSigScoreMatrix <- matrix(unlist(bplapply(randomSigData, function(rsig) { rank(rsig@scores, ties.method="average") },
  										BPPARAM=BPPARAM)), nrow=N_SAMPLES, ncol=length(randomSigData))


  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "rSigScoreMatrix")


  ### build one hot matrix for factors
  factorSigs <- list()
  for (s in precomputedFactor) {
    fValues <- s@scores
    fLevels <- unique(fValues)
    factorFreq <- matrix(0L, ncol=length(fLevels))
    factorMatrix <- matrix(0L, nrow=N_SAMPLES, ncol=length(fLevels))

    factList <- lapply(fLevels, function(fval) {
				factorMatrixRow <- matrix(0L, nrow=N_SAMPLES, ncol=1)
				equal_ii <- which(fValues == fval)
				factorMatrixRow[equal_ii] <- 1
				return(list(length(equal_ii) / length(fValues), factorMatrixRow))
			})
	factorFreq <- lapply(factList, function(x) return(x[[1]]))
	factorMatrix <- matrix(unlist(lapply(factList, function(x) return(x[[2]]))), nrow=N_SAMPLES, ncol=length(fLevels))

    factorSigs[[s@name]] <- list(fLevels, factorFreq, factorMatrix)
  }

  timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  tRows <- c(tRows, "Factor Matrix")

  message("Evaluating signatures against projections...")

  i <- 1
  projnames <- names(projections)
  for (proj in projections) {
  	weights <- computeKNNWeights(proj, K=round(sqrt(NCOL(proj@pData))), BPPARAM)
  	neighborhoodPrediction <- Matrix::crossprod(weights, sigScoreMatrix)

  	## Neighborhood dissimilatory score = |actual - predicted|
  	dissimilarity <- abs(sigScoreMatrix - neighborhoodPrediction)
  	medDissimilarity <- as.matrix(apply(dissimilarity, 2, median))

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "Dissimilarity"))


    # Calculate scores for random signatures
    randomNeighborhoodPrediction <- Matrix::crossprod(weights, randomSigScoreMatrix)
    randomDissimilarity <- abs(randomSigScoreMatrix - randomNeighborhoodPrediction)
    randomMedDissimilarity <- as.matrix(apply(randomDissimilarity, 2, median))

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "RandDissimilarity"))

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

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "GroupNumGenes"))

      bgStat <- matrix(unlist(bplapply(names(backgrounds), function(x) {
  			mu_x <- mean(backgrounds[[x]])
  			std_x <- biasedVectorSD(as.matrix(backgrounds[[x]]))
  			return(list(as.numeric(x), mu_x, std_x))
  		  }, BPPARAM=BPPARAM)), nrow=length(names(backgrounds)), ncol=3, byrow=T)

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "bgStat"))


  	mu <- matrix(unlist(bplapply(spRows, function(x) {
  			numG <- x@numGenes
  			row_i <- which.min(abs(numG - bgStat[,1]))
  			return(bgStat[row_i, 2])
  		  }, BPPARAM=BPPARAM)), nrow=nrow(medDissimilarity), ncol=ncol(medDissimilarity))

  	sigma <- matrix(unlist(bplapply(spRows, function(x) {
  			numG <- x@numGenes
  			row_i <- which.min(abs(numG - bgStat[,1]))
  			return(bgStat[row_i,3])
  		  }, BPPARAM=BPPARAM)), nrow=nrow(medDissimilarity), ncol=ncol(medDissimilarity))

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "Mu/Sigma"))


      #Create CDF function for medDissmilarityPrime and apply CDF function to medDissimilarityPrime pointwise
      pValues <- pnorm( ((medDissimilarity - mu) / sigma))

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "PValues"))

      sigProjMatrix[,i] <- 1 - (medDissimilarity / N_SAMPLES)
  	sigProjMatrix_P[,i] <- pValues

      # Calculate significance for precomputed numerical signatures
      # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)
      #j <- 1
      #for (s in precomputedNumerical) {

      #  sigScores <- rank(s@scores)
      #  if (all(sigScores == sigScores[1])) {
      #    pnumSigProjMatrix[j,i] <- 0.0
      #    pnumSigProjMatrix_P[j,i] <- 1.0
      #    next
      #  }
      #  sigPredictions <- crossprod(weights, sigScores)
      #  r <- roc.area(sigScores, sigPredictions)
      #  a <- r$A
        #dissimilarity <- abs(sigScores - sigPredictions)
        #medDissimilarity <- as.matrix(apply(dissimilarity, 2, median))

        #Compute a background for numerical signatures
        #NUM_REPLICATES <- 10000
        #randomSigValues <- get_bg_dist(N_SAMPLES, NUM_REPLICATES)
        #bgValues <- as.vector(t(sigScores))[randomSigValues]
        #randomPredictions <- (weights %*% bgValues)
        #randR <- roc(sigScores, sigPredictions)
        #randA <- randR$auc

        #rDissimilarity <- abs(bgValues <- randomPredictions)
        #randomScores <- as.matrix(apply(rDissimilarity, 2, median))

        #mu <- mean(randomScores)
        #sigma <- biasedVectorSD(random_scores)
        #if (sigma != 0) {
        #  p_value <- pnorm((medDissimilarity - mu) / sigma)
        #} else {
        #  p_value <- 1.0
        #}

        #pnumSigProjMatrix[j,i] <- 1 - (medDissimilarity / N_SAMPLES)
        #pnumSigProjMatrix_P[j,i] <- p_value


      #  pnumSigProjMatrix[j,i] <- a
      #  pnumSigProjMatrix_P[j,i] <- r$p.value

      #  j <- j + 1
      #}


      ## Calculate signficance for Factor signatures

  	factorSigProjList <- lapply(factorSigs, function(x) {
  			fLevels <- x[[1]]
  			factorFreq <- x[[2]]
  			fMatrix <- x[[3]]

  			if (1 %in% factorFreq) {
  				return(list(0.0, 1.0))
  			}

  			N_LEVELS <- length(fLevels)
  			factorPredictions <- Matrix::crossprod(weights, fMatrix)

  			labels <- apply(fMatrix, 1, which.max)

  			# Get predicted index and corresponding probability
  			preds_ii <- apply(factorPredictions, 1, which.max)
  			preds <- factorPredictions[cbind(1:nrow(factorPredictions), preds_ii)]

  			r <- multiclass.roc(labels, preds_ii, levels=seq(1, length(fLevels)))
  			a <- r$auc

  			# Calculate the p value for precomputed signature
  			krList <- list()
  			for (k in 1:length(fLevels)) {
  				krList <- c(krList, list(preds_ii[labels==k]))
  			}

  			krTest <- kruskal.test(krList)

  			return(list(a, krTest$p.value))
  		  })

  	factorSigProjMatrix[,i] <- unlist(lapply(factorSigProjList, function(x) x[[1]]))
  	factorSigProjMatrix_P[,i] <- unlist(lapply(factorSigProjList, function(x) x[[2]]))

  	timingList <- rbind(timingList, c(difftime(Sys.time(), ptm, units="secs")))
  	tRows <- c(tRows, paste0(proj@name, "FactorSigs"))

    i <- i+1
  }

  # Concatenate Factor Sig-proj entries back in
  sigProjMatrix <- rbind(sigProjMatrix, factorSigProjMatrix, pnumSigProjMatrix)
  sigProjMatrix_P <- rbind(sigProjMatrix_P, factorSigProjMatrix_P, pnumSigProjMatrix_P)
  spRowLabels <- c(spRowLabels, spRowLabelsFactors, spRowLabelsPNum)

  original_shape <- dim(sigProjMatrix_P)
  sigProjMatrix_P <- matrix(p.adjust(sigProjMatrix_P, method="BH"), nrow=nrow(sigProjMatrix_P), ncol=ncol(sigProjMatrix_P))
  sigProjMatrix_P[sigProjMatrix_P == 0] <- 10^(-300)

  sigProjMatrix_P <- as.matrix(log10(sigProjMatrix_P))

  colnames(sigProjMatrix_P) <- projnames
  rownames(sigProjMatrix_P) <- spRowLabels

  colnames(sigProjMatrix) <- projnames
  rownames(sigProjMatrix) <- spRowLabels

  timingList <- as.matrix(timingList)
  rownames(timingList) <- tRows

  #return(list(list(spRowLabels, spColLabels, sigProjMatrix, sigProjMatrix_P), timingList))
  return(list(sigNames = spRowLabels, projNames = spColLabels,
              sigProjMatrix = sigProjMatrix, pVals = sigProjMatrix_P))
}

#' Clusters signatures according to the rank sum
#'
#' @param sigList List of signatures
#' @param sigMatrix Matrix of signatures scores, NUM_SIGNATURES x NUM_SAMPLES
#' @param k Number of clusters to generate
#' @return List of clusters involving computed signatures
#' @return List of clusters involving precomptued signatures
clusterSignatures <- function(sigList, sigMatrix, pvals, k=10) {

  precomputed <- lapply(sigList, function(x) x@isPrecomputed)
  insignificant <- apply(pvals, 1, function(x) min(x) < -1.3)

  keep = names(which(insignificant == T))

  # Cluster computed signatures and precomputed signatures separately
  computedSigsToCluster <- names(precomputed[which(precomputed==F)])
  computedSigMatrix <- sigMatrix[computedSigsToCluster,,drop=F]

  computedSigMatrix <- computedSigMatrix[keep,,drop=F]

  compcls <- list()
  maxcls <- 1
  if (nrow(computedSigMatrix) > 1) {
  	r <- as.matrix(t(apply(computedSigMatrix, 1, function(x) rank(x, ties.method="average"))))

	compkm <- densityMclust(r)
	compcls <- as.list(compkm$classification)
	compcls <- compcls[order(unlist(compcls), decreasing=F)]

	maxcls <- max(unlist(compcls))
  }

  compcls[names(which(insignificant == F))] <- maxcls + 1

  # Don't actually cluster Precomputed Signatures -- just return in a list.
  precompcls <- list()
  if (length(which(precomputed==T)) > 0) {
  	  precomputedSigsToCluster <- names(precomputed[which(precomputed==T)])
  	  precompcls[precomputedSigsToCluster] <- 1
  }


  output <- list(compcls, precompcls)
  names(output) <- c("Computed", "Precomputed")

  return(output)
}

calcMutualInformation <- function(ss, pc) {

	minVal <- min(min(ss), min(pc))
	maxVal <- max(max(ss), max(pc))

	f1 <- hist(ss, seq(minVal-0.1, maxVal+0.1, by=0.1), plot=F)$counts
	f2 <- hist(pc, seq(minVal-0.1, maxVal+0.1, by=0.1), plot=F)$counts

	jdist <- cbind(f1, f2)

	h1 <- entropy.plugin(rowSums(jdist))
	h2 <- entropy.plugin(colSums(jdist))
	h12 <- entropy.plugin(jdist)
	mi <- h1 + h2 - h12

	return(mi)
}

calcPearsonCorrelation <- function(ss, pc) {

	pc <- rcorr(ss, pc, type="pearson")

	return(pc[[1]][1,2])

}

