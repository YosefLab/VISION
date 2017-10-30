#' Initialize a new Signature object.
#' Should not be called directly, instead use the `new` syntax
#'
#' @param sigDict Named list of signs for each gene in the signature
#' @param name Name of the signature
#' @param source File from which this signature was read from
#' @param metaData Metadata pertinent to signature
#' @param isPrecomputed If TRUE indicates that this signature was precomputed.
#' Else not precomputed. Default is FALSE.
#' @param isFactor If TRUE indicates that this signature is a Factor, else not
#' a factor. Default is FALSE.
#' @param cluster Number representing which cluster this signature is a part of.
#' Default is 0.
#' @return Signature object
Signature <- function(sigDict, name, source, metaData="",
                      isPrecomputed=FALSE,
                      isFactor=FALSE, cluster=0) {
  if (missing(sigDict)) {
    stop("Missing sigDict information.")
  } else if (missing(name)) {
    stop("Missing signature name.")
  } else if (missing(source)) {
    stop("Missing source file.")
  }

  .Object <- new("Signature", sigDict=sigDict, name=name, source=source,
                 metaData=metaData, isPrecomputed=isPrecomputed,
                 isFactor=isFactor, cluster=cluster)

  return(.Object)
}

#' Create a user-defined gene signature
#'
#' @param name the name of the signature
#' @param sigData a named vector where the names correspond to genes in the
#' data and the values are either `1.0` for up-regulated (or positive) genes,
#' and `-1.0` for down regulated (negtive) genes.
#' @param metadata metadata that is relevent to the signature. [Default:NULL]
#' @export
#' @return a Signature object
#' @examples
#' sig <- createGeneSignature(name="signature", c(gene1=1.0,gene17=1.0,
#'                                                      gene4=1.0, gene31=-1.0,
#'                                                      gene3=-1.0, gene9=1.0))
#' sig2 <- createGeneSignature(name="signature", c(gene18=1.0,gene29=1.0,
#'                                                      gene400=-1.0, gene1=-1.0,
#'                                                      gene7=1.0, gene9=1.0))
createGeneSignature <- function(name, sigData, metadata="") {
  return(new("Signature", sigDict=sigData, name=name, metaData=metadata,
             source="user-defined"))
}

#' Method for testing whether or not signatures are equal to one another
#' (if all signs, genes and name are the same.)
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


#' calculate background signature scores
#'
#' For each signature-cell pair, compute a score that captures the level of
#' correspondence between the cell and the signature.
#' To estimate significance of these scores, a set of random gene signatures is
#' generated to create a null distribution
#'
#' @param object the FastProject object
#' @param num the number of signatures to generate
#' @param BPPARAM the parallelization backend to use for these computations
#' @return numeric matrix of random signature scores SIGNATURES x CELLS
calculateSignatureBackground <- function(object, num, BPPARAM=NULL) {

    if(is.null(BPPARAM)) BPPARAM <- SerialParam()

    message("Computing background distribution for signature scores...")
    # Construct random signatures for background distribution
    sigSizes <- lapply(object@sigData, function(s) length(s@sigDict))
    randomSigs <- generatePermutationNull(num, object@exprData, sigSizes)

    normExpr <- getNormalizedCopy(object@exprData, object@sig_norm_method)

    ## Compute signature scores for random signatures generated
    randomSigScores <- bplapply(randomSigs, function(s) {
        singleSigEval(s, object@sig_score_method, normExpr,
                      object@weights,
                      object@min_signature_genes)
    } ,BPPARAM=BPPARAM)
    names(randomSigScores) <- names(randomSigs)

    ## Remove random signatures that didn't compute correctly
    toRemove <- vapply(randomSigScores, is.null, TRUE)
    randomSigScores <- randomSigScores[!toRemove]

    return(randomSigScores)

}

#' Generate random signatures for a null distribution by permuting the data
#' @param num the number of signatures to generate
#' @param eData the data to use for the permutations
#' @param sigSizes the sizes of the true signatures, to inform selection of the
#' random signature sizes
#' @return a vector of random Signature objects
generatePermutationNull <- function(num, eData, sigSizes) {
  randomSigs <- c()
  randomSizes <- unlist(findRepSubset(sigSizes))
  for (size in randomSizes) {
    for (j in 1:num) {
      newSigGenes <- sample(rownames(getExprData(eData)), min(size, nrow(getExprData(eData))))
      newSigSigns <- rep(1, size)
      names(newSigSigns) <- newSigGenes
      newSig <- Signature(newSigSigns, paste0("RANDOM_BG_", size, "_", j), 'x')
      randomSigs <- c(randomSigs, newSig)
    }
  }
  return(randomSigs)
}

BG_DIST <- matrix(0L, nrow=0, ncol=0)

#' Generates, or if already generated retrieves, a random Background Distribution
#' @importFrom stats rnorm
#' @param N_SAMPLES Number of samples to generate the distribution for
#' @param NUM_REPLICATES Number of replicates to generate for the background distribution
#' @return Random matrix with dimensions N_SAMPLES x NUM_REPLICATES with row ordered in ascending order
getBGDist <- function(N_SAMPLES, NUM_REPLICATES) {

  if (nrow(BG_DIST) != N_SAMPLES || ncol(BG_DIST) != NUM_REPLICATES) {
    BG_DIST <- matrix(rnorm(N_SAMPLES*NUM_REPLICATES),
                      nrow=N_SAMPLES, ncol=NUM_REPLICATES)
    BG_DIST <- apply(BG_DIST, 2, order)
  }

  return(BG_DIST)

}

#' Evaluates the significance of each signature vs. each projection.
#' @importFrom entropy entropy.plugin
#' @importFrom stats p.adjust
#' @param projections Maps projections to their spatial coordinates for each
#' sample
#' @param sigScoresData List of SignatureScores Object, mapping signature
#' names to their value at each coordinate
#' @param randomSigData List of SignatureScores Object, mapping randomly
#' generated signatures to scores to be compared with the real signature scores.
#' @param BPPARAM the parallelization backend to use
#' @return list:
#' \itemize{
#'     \item sigProbMatrix: the matrix of signature-projection consistency scores
#'     \item pVals: pvalues for the scores
#' }
sigsVsProjections <- function(projections, sigScoresData,
                              randomSigData, BPPARAM=bpparam()) {

  N_SAMPLES <- length(sigScoresData[[1]]@sample_labels)

  # Build a matrix of all non-precomputed signatures
  geneSigs <- sigScoresData[vapply(sigScoresData,
                                   function(x) !x@isPrecomputed,
                                   TRUE)]

  sigScoreMatrix <- do.call(
               cbind,
               lapply(geneSigs, function(sig) sig@scores)
  )
  sampleLabels <- rownames(sigScoreMatrix)
  sigScoreMatrix <- matrixStats::colRanks(
             sigScoreMatrix, preserveShape=TRUE, ties.method="average"
  )
  colnames(sigScoreMatrix) <- names(geneSigs)
  rownames(sigScoreMatrix) <- sampleLabels


  # Build a matrix of random background signatures
  randomSigScoreMatrix <- do.call(
               cbind,
               lapply(randomSigData, function(rsig) rsig@scores)
  )
  sampleLabels <- rownames(randomSigScoreMatrix)
  randomSigScoreMatrix <- matrixStats::colRanks(
             randomSigScoreMatrix, preserveShape=TRUE, ties.method="average"
  )

  # TODO: Set names of randomSigData in analyze and use
  colnames(randomSigScoreMatrix) <- vapply(randomSigData,
                                           function(x) x@name, "")
  rownames(randomSigScoreMatrix) <- sampleLabels


  message("Evaluating signatures against projections...")

  sigProjMatrix <- data.frame()
  sigProjMatrix_P <- data.frame()

  for (proj in projections) {
    weights <- computeKNNWeights(proj, K=round(sqrt(NCOL(proj@pData))), BPPARAM)

    svp_n <- sigsVsProjection_n(geneSigs, sigScoreMatrix,
                                randomSigData, randomSigScoreMatrix, weights)
    svp_pcn <- sigsVsProjection_pcn(sigScoresData, weights)
    svp_pcf <- sigsVsProjection_pcf(sigScoresData, weights)

    consistency = c(svp_n$consistency,
                    svp_pcn$consistency,
                    svp_pcf$consistency)

    consistency <- as.data.frame(consistency)
    colnames(consistency) <- c(proj@name)

    pvals = c(svp_n$pvals,
              svp_pcn$pvals,
              svp_pcf$pvals)

    pvals <- as.data.frame(pvals)
    colnames(pvals) <- c(proj@name)

    sigProjMatrix <- merge(sigProjMatrix, consistency,
                           by='row.names', all=TRUE)

    rownames(sigProjMatrix) = sigProjMatrix$Row.names
    sigProjMatrix$Row.names <- NULL

    sigProjMatrix_P <- merge(sigProjMatrix_P, pvals,
                             by='row.names', all=TRUE)

    rownames(sigProjMatrix_P) = sigProjMatrix_P$Row.names
    sigProjMatrix_P$Row.names <- NULL

  }

  # FDR-correct

  sigProjMatrix_Padj <- matrix(
      p.adjust(
          as.matrix(sigProjMatrix_P),
          method="BH"
      ),
    nrow=nrow(sigProjMatrix_P),
    ncol=ncol(sigProjMatrix_P)
  )
  colnames(sigProjMatrix_Padj) <- colnames(sigProjMatrix_P)
  rownames(sigProjMatrix_Padj) <- rownames(sigProjMatrix_P)

  sigProjMatrix_Padj[sigProjMatrix_Padj == 0] <- 10^(-300)

  sigProjMatrix <- as.matrix(sigProjMatrix)
  sigProjMatrix_Padj <- as.matrix(log10(sigProjMatrix_Padj))

  return(list(sigProjMatrix = sigProjMatrix, pVals = sigProjMatrix_Padj))
}

#' Evaluates the significance of each numeric signature vs. a
#' single projections weights
#' @param sigData a list of SignatureScores objects representing input user
#' signatures
#' @param sigScoreMatrix a numeric matrix (N_SAMPLES x N_SIGNATURES) containing
#' the signature scores for every (signature, sample)
#' @param randomSigData a list of SignatureScores objects representing a suitable
#' background distribution to go with `sigData`
#' @param randomSigScoreMatrix a numeric matrix (N_SAMPLES x N_SIGNATURES) containing
#' the signature scores associated wiht signatures in randomSigData
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#' }
sigsVsProjection_n <- function(sigData, sigScoreMatrix,
                               randomSigData, randomSigScoreMatrix, weights){
  # Calculate significance for precomputed numerical signatures
  # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)

  N_SAMPLES = nrow(weights)

  neighborhoodPrediction <- weights %*% sigScoreMatrix

  ## Neighborhood dissimilatory score = |actual - predicted|
  dissimilarity <- abs(sigScoreMatrix - neighborhoodPrediction)
  medDissimilarity <- as.matrix(apply(dissimilarity, 2, median))

  # Calculate scores for random signatures
  randomNeighborhoodPrediction <- weights %*% randomSigScoreMatrix
  randomDissimilarity <- abs(randomSigScoreMatrix -
                               randomNeighborhoodPrediction)
  randomMedDissimilarity <- as.matrix(apply(randomDissimilarity, 2,
                                            median))

  # Group by number of genes
  backgrounds <- list()
  k <- 1
  for (rsig in randomSigData) {
    numGenes <- as.character(rsig@numGenes)

    if (!(numGenes %in% names(backgrounds))) {
      backgrounds[[numGenes]] <- c(randomMedDissimilarity[k])
    } else {
      backgrounds[[numGenes]] <- c(backgrounds[[numGenes]],
                                   randomMedDissimilarity[k])
    }
    k <- k + 1
  }

  # TODO: Do the lookup here using the column names of randomSigScoreMatrix
  bgStat <- matrix(unlist(lapply(names(backgrounds), function(x) {
    mu_x <- mean(backgrounds[[x]])
    std_x <- biasedVectorSD(as.matrix(backgrounds[[x]]))
    return(list(as.numeric(x), mu_x, std_x))
  })), nrow=length(names(backgrounds)),
  ncol=3, byrow=TRUE)


  mu <- matrix(unlist(lapply(colnames(sigScoreMatrix), function(x) {
    numG <- sigData[[x]]@numGenes
    row_i <- which.min(abs(numG - bgStat[,1]))
    return(bgStat[row_i, 2])
  })), nrow=nrow(medDissimilarity),
  ncol=ncol(medDissimilarity))

  sigma <- matrix(unlist(lapply(colnames(sigScoreMatrix), function(x) {
    numG <- sigData[[x]]@numGenes
    row_i <- which.min(abs(numG - bgStat[,1]))
    return(bgStat[row_i,3])
  })), nrow=nrow(medDissimilarity),
  ncol=ncol(medDissimilarity))


  #Create CDF function for medDissmilarityPrime and apply CDF function to
  #medDissimilarityPrime pointwise
  pvals <- pnorm( ((medDissimilarity - mu) / sigma))

  consistency <- 1 - (medDissimilarity / N_SAMPLES)

  #Reduce from matrix col to char-vector
  pvals <- pvals[,1]
  consistency <- consistency[,1]


  return(list(consistency=consistency, pvals=pvals))
}

#' Evaluates the significance of each precomputed numeric signature vs. a
#' single projections weights
#' @importFrom stats pnorm
#' @param sigScoresData List of SignatureScores Object, mapping signature
#' names to their value at each coordinate
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#' }
sigsVsProjection_pcn <- function(sigScoresData, weights){
  # Calculate significance for precomputed numerical signatures
  # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)

  N_SAMPLES = nrow(weights)

  consistency <- numeric()
  pvals <- numeric()

  for (s in sigScoresData) {

    if( !s@isPrecomputed || s@isFactor) {
      next
    }

    sigScores <- rank(s@scores, ties.method="average")
    if (all(sigScores == sigScores[1])) {
      consistency[s@name] <- 0
      pvals[s@name] <- 1.0
      next
    }

    sigPredictions <- weights %*% sigScores
    dissimilarity <- abs(sigScores - sigPredictions)
    medDissimilarity <- median(dissimilarity)

    #Compute a background for numerical signatures
    NUM_REPLICATES <- 3000
    bgValues <- replicate(NUM_REPLICATES, sample(sigScores))
    randomPredictions <- weights %*% bgValues

    rDissimilarity <- abs(bgValues - randomPredictions)
    randomScores <- as.matrix(apply(rDissimilarity, 2, median))


    mu <- mean(randomScores)
    sigma <- biasedVectorSD(randomScores)

    if (sigma != 0) {
      p_value <- pnorm((medDissimilarity - mu) / sigma)
    } else {
      p_value <- 1.0
    }

    consistency[s@name] <- 1 - (medDissimilarity / N_SAMPLES)
    pvals[s@name] <- p_value

  }

  return(list(consistency=consistency, pvals=pvals))
}

#' Evaluates the significance of each precomputed factor signature vs. a
#' single projections weights
#' @importFrom stats kruskal.test
#' @param sigScoresData List of SignatureScores Object, mapping signature
#' names to their value at each coordinate
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#' }
sigsVsProjection_pcf <- function(sigScoresData, weights){
  # Calculate significance for precomputed numerical signatures
  # This is done separately because there are likely to be many repeats (e.g. for a time coordinate)

  N_SAMPLES = nrow(weights)

  consistency <- numeric()
  pvals <- numeric()

  for (s in sigScoresData) {

    if( !s@isPrecomputed || !s@isFactor) {
      next
    }

    ### build one hot matrix for factors

    fValues <- s@scores
    fLevels <- levels(fValues)
    factorFreq <- matrix(0L, ncol=length(fLevels))
    factorMatrix <- matrix(0L, nrow=N_SAMPLES, ncol=length(fLevels))

    factList <- lapply(fLevels, function(fval) {
      factorMatrixRow <- matrix(0L, nrow=N_SAMPLES, ncol=1)
      equal_ii <- fValues == fval
      factorMatrixRow[equal_ii] <- 1
      return(list(sum(equal_ii) / length(fValues), factorMatrixRow))
    })
    factorFreq <- lapply(factList, function(x) return(x[[1]]))
    factorMatrix <- matrix(unlist(lapply(factList, function(x) return(x[[2]]))),
                           nrow=N_SAMPLES, ncol=length(fLevels))

    if (1 %in% factorFreq) {
      consistency[s@name] <- 0.0
      pvals[s@name] <- 1.0
      next
    }

    factorPredictions <- weights %*% factorMatrix

    labels <- apply(factorMatrix, 1, which.max)

    # Get predicted index and corresponding probability
    preds_ii <- apply(factorPredictions, 1, which.max)
    preds <- factorPredictions[cbind(1:nrow(factorPredictions), preds_ii)]

    # Calculate the p value for precomputed signature
    krList <- list()
    for (k in 1:length(fLevels)) {
      krList <- c(krList, list(preds_ii[labels==k]))
    }

    krTest <- kruskal.test(krList)

    consistency[s@name] <- krTest$statistic
    pvals[s@name] <- krTest$p.value

  }

  return(list(consistency=consistency, pvals=pvals))
}

#' Clusters signatures according to the rank sum
#' @importFrom mclust Mclust
#' @param sigList List of signatures
#' @param sigMatrix data.frame of signatures scores, NUM_SAMPLES x NUM_SIGNATURES
#' @param pvals the corresponding P-values for each score,
#' NUM_SIGNATURES x NUM_SAMPLES
#' @return a list:
#' \itemize{
#'     \item Computed: a list of clusters of computed signatures
#'     \item Precomputed: a list of clusters of precomputed signatures
#' }
clusterSignatures <- function(sigList, sigMatrix, pvals) {

  precomputed <- vapply(sigList, function(x) x@isPrecomputed, TRUE)
  significant <- apply(pvals, 1, function(x) min(x) < -1.3)

  comp_names <- names(precomputed)[!precomputed]
  signif_computed <- significant[comp_names]

  keep_computed = names(signif_computed)[signif_computed]

  # Cluster computed signatures and precomputed signatures separately
  computedSigMatrix <- sigMatrix[, keep_computed,drop=FALSE]

  compcls <- list()
  maxcls <- 1
  if (ncol(computedSigMatrix) > 1) {

    r <- colRanks(
            as.matrix(computedSigMatrix),
            ties.method="average",
            preserveShape=TRUE
        )

    compkm <- Mclust(t(r))

    compcls <- as.list(compkm$classification)
    names(compcls) = colnames(computedSigMatrix)
    compcls <- compcls[order(unlist(compcls), decreasing=FALSE)]

    maxcls <- max(unlist(compcls))
  }

  compcls[names(signif_computed)[!signif_computed]] <- maxcls + 1

  # Don't actually cluster Precomputed Signatures -- just return in a list.
  precompcls <- list()
  if (any(precomputed)) {
    precomputedSigsToCluster <- names(precomputed)[precomputed]
    precompcls[precomputedSigsToCluster] <- 1
  }

  return(list(Computed = compcls, Precomputed = precompcls))
}

#' Compute pearson correlation between signature scores and principle components
#' @importFrom Hmisc rcorr
#' @param ss signature scores
#' @param pc principle components
#' @return pearson correlaton coeffcients
calcPearsonCorrelation <- function(ss, pc) {

  pc <- rcorr(ss, pc, type="pearson")

  return(pc[[1]][1,2])

}
