#' Initialize a new Signature object.
#' Should not be called directly, instead use the `new` syntax
#'
#' @param sigDict Named list of signs for each gene in the signature
#' @param name Name of the signature
#' @param source File from which this signature was read from
#' @param metaData Metadata pertinent to signature
#' @param isMeta If TRUE indicates that this signature was derived from
#' meta data. Default is FALSE.
#' @param isFactor If TRUE indicates that this signature is a Factor, else not
#' a factor. Default is FALSE.
#' @return Signature object
Signature <- function(sigDict, name, source, metaData="",
                      isMeta=FALSE,
                      isFactor=FALSE) {
  if (missing(sigDict)) {
    stop("Missing sigDict information.")
  } else if (missing(name)) {
    stop("Missing signature name.")
  } else if (missing(source)) {
    stop("Missing source file.")
  }

  names(sigDict) <- vapply(
                       names(sigDict), toupper,
                       "", USE.NAMES = FALSE
                       )

  .Object <- new("Signature", sigDict=sigDict, name=name, source=source,
                 metaData=metaData, isMeta=isMeta,
                 isFactor=isFactor)

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
#' @return A list with two items:
#'
#'   randomSigs: a list of lists of SignatureScore objects.  Each sub-list
#'     represents permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
calculateSignatureBackground <- function(object, num) {

    # Construct random signatures for background distribution
    out <- generatePermutationNull(num, object@exprData, object@sigData)
    randomSigs <- out$randomSigs
    sigAssignments <- out$sigAssignments

    normExpr <- getNormalizedCopy(object@exprData, object@sig_norm_method)

    randomSigScores <- batchSigEval(
                            unlist(randomSigs, recursive = FALSE),
                            object@sig_score_method,
                            normExpr, object@weights,
                            object@min_signature_genes)

    # randomSigScores is list of SignatureScore
    # need to make a list of lists of SignatureScores, grouping
    # by background group number
    randomSigScoresGroups <- list()
    for (groupName in names(randomSigs)) {
        sigsInGroup <- names(randomSigs[[groupName]])
        randomSigScoresGroups[[groupName]] <- randomSigScores[sigsInGroup]
    }

    return(list(
                randomSigs = randomSigScoresGroups,
                sigAssignments = sigAssignments))

}

#' Generate random signatures for a null distribution by permuting the data
#'
#' @param num the number of signatures to generate
#' @param eData the data to use for the permutations
#' @param sigData list of signature objects
#' random signature sizes
#' @return A list with two items:
#'
#'   randomSigs: a list of lists of Signature objects.  Each sub-list represents
#'     permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
generatePermutationNull <- function(num, eData, sigData) {

  exp_genes <- rownames(eData)

  # Remove meta
  meta <- vapply(sigData, function(s) s@isMeta, TRUE)
  sigData <- sigData[!meta]

  sigSize <- vapply(sigData, function(s) {
                        sum(names(s@sigDict) %in% exp_genes)
                    }
                 , 1)

  sigSize <- log10(sigSize)
  sigBalance <- vapply(sigData, function(s) {
                           positive <- sum(s@sigDict >= 0)
                           balance <- positive / length(s@sigDict)
                           return(balance)
    }, 1)

  sigBalance[sigBalance < 0.5] <- 1 - sigBalance[sigBalance < 0.5]

  #sigVars = cbind(sigSize, sigBalance)
  sigVars <- cbind(sigSize) # TODO: incorporate the sigBalance here

  n_components <- 5 # TODO: choose number of components better
  if (nrow(unique(sigVars)) < n_components){
      n_components <- nrow(unique(sigVars))
  }

  km <- kmeans(sigVars, n_components)
  centers <- km$centers
  clusters <- as.factor(km$cluster)

  message("Creating ", nrow(centers),
          " background signature groups with the following parameters:")
  centers[, "sigSize"] <- round(10 ** centers[, "sigSize"])
  print(centers) # How do I do this with 'message'??
  message("  signatures per group: ", num)

  randomSigs <- list()

  for (cluster_i in rownames(centers)) {

    clusterSigs <- list()

    size <- centers[cluster_i, "sigSize"]
    #balance <- centers[cluster_i,"sigBalance"]
    balance <- 1

    for (j in 1:num) {
      newSigGenes <- sample(exp_genes, min(size, length(exp_genes)))
      newSigSigns <- sample(c(1, -1), length(newSigGenes),
                            replace = TRUE, prob = c(balance, 1 - balance))
      names(newSigSigns) <- newSigGenes
      newSig <- Signature(newSigSigns, paste0("RANDOM_BG_", size, "_", j), "x")
      clusterSigs[[newSig@name]] <- newSig
    }

    randomSigs[[cluster_i]] <- clusterSigs
  }

  return(list(randomSigs = randomSigs, sigAssignments = clusters))
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

#' Evaluates the significance of each signature in each cluster
#'
#' @importFrom stats p.adjust
#' @param latentSpace numeric matrix N_Cells x N_Components
#' @param sigScoresData List of SignatureScores Object, mapping signature
#' names to their value at each coordinate
#' @param metaData data.frame of meta-data for cells
#' @param randomSigData A list with two items:
#'
#'   randomSigs: a list of lists of SignatureScore objects.  Each sub-list
#'     represents permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
#' @param clusters a factor vector denoting which cluster each cell belongs to
#' @return list:
#' \itemize{
#'     \item sigProbMatrix: the matrix of signature-projection consistency scores
#'     \item pVals: pvalues for the scores
#' }
sigsVsProjectionsClusters <- function(latentSpace, sigScoresData, metaData,
                              randomSigData, clusters) {

  sigProjMatrix <- data.frame()
  sigProjMatrix_P <- data.frame()
  sigProjMatrix_Pemp <- data.frame()

  proj <- Projection("All", pData = t(latentSpace))

  # First evaluate things for all clusters in the latent space
  weights <- computeKNNWeights(proj, K = round(sqrt(NCOL(proj@pData))))

  svp_n <- sigsVsProjection_n(sigScoresData,
                              randomSigData, weights)
  svp_pcn <- sigsVsProjection_pcn(metaData, weights)
  svp_pcf <- sigsVsProjection_pcf(metaData, weights)

  consistency <- c(svp_n$consistency,
                  svp_pcn$consistency,
                  svp_pcf$consistency)

  consistency <- as.data.frame(consistency)
  colnames(consistency) <- c(proj@name)

  pvals <- c(svp_n$pvals,
            svp_pcn$pvals,
            svp_pcf$pvals)

  emp_pvals <- c(svp_n$empvals,
                svp_pcn$pvals,
                svp_pcf$pvals)

  pvals <- as.data.frame(pvals)
  colnames(pvals) <- c(proj@name)

  emp_pvals <- as.data.frame(emp_pvals)
  colnames(emp_pvals) <- c(proj@name)

  sigProjMatrix <- merge(sigProjMatrix, consistency,
                         by = "row.names", all = TRUE)

  rownames(sigProjMatrix) <- sigProjMatrix$Row.names
  sigProjMatrix$Row.names <- NULL

  sigProjMatrix_P <- merge(sigProjMatrix_P, pvals,
                           by = "row.names", all = TRUE)

  sigProjMatrix_Pemp <- merge(sigProjMatrix_Pemp, emp_pvals,
                              by = "row.names", all = TRUE)

  rownames(sigProjMatrix_P) <- sigProjMatrix_P$Row.names
  rownames(sigProjMatrix_Pemp) <- sigProjMatrix_Pemp$Row.names
  sigProjMatrix_P$Row.names <- NULL
  sigProjMatrix_Pemp$Row.names <- NULL

  for (cluster in unique(clusters)) {
    cluster_cells <- names(clusters)[clusters == cluster]
    proj <- Projection(as.character(cluster),
                       pData = t(latentSpace[cluster_cells, ]))


    weights <- computeKNNWeights(proj, K=round(sqrt(NCOL(proj@pData))))

    svp_n <- sigsVsProjection_n(sigScoresData, randomSigData,
                                weights, cells = cluster_cells)
    svp_pcn <- sigsVsProjection_pcn(metaData, weights, cells = cluster_cells)
    svp_pcf <- sigsVsProjection_pcf(metaData, weights, cells = cluster_cells)

    consistency <- c(svp_n$consistency,
                    svp_pcn$consistency,
                    svp_pcf$consistency)

    consistency <- as.data.frame(consistency)
    colnames(consistency) <- c(proj@name)

    pvals <- c(svp_n$pvals,
              svp_pcn$pvals,
              svp_pcf$pvals)

    emp_pvals <- c(svp_n$empvals,
                  svp_pcn$pvals,
                  svp_pcf$pvals)

    pvals <- as.data.frame(pvals)
    colnames(pvals) <- c(proj@name)

    emp_pvals <- as.data.frame(emp_pvals)
    colnames(emp_pvals) <- c(proj@name)

    sigProjMatrix <- merge(sigProjMatrix, consistency,
                           by = "row.names", all = TRUE)

    rownames(sigProjMatrix) <- sigProjMatrix$Row.names
    sigProjMatrix$Row.names <- NULL

    sigProjMatrix_P <- merge(sigProjMatrix_P, pvals,
                             by = "row.names", all = TRUE)

    sigProjMatrix_Pemp <- merge(sigProjMatrix_Pemp, emp_pvals,
                                by = "row.names", all = TRUE)

    rownames(sigProjMatrix_P) <- sigProjMatrix_P$Row.names
    rownames(sigProjMatrix_Pemp) <- sigProjMatrix_Pemp$Row.names
    sigProjMatrix_P$Row.names <- NULL
    sigProjMatrix_Pemp$Row.names <- NULL

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

  sigProjMatrix_Pempadj <- matrix(
      p.adjust(
          as.matrix(sigProjMatrix_Pemp),
          method="BH"
      ),
    nrow=nrow(sigProjMatrix_Pemp),
    ncol=ncol(sigProjMatrix_Pemp)
  )
  colnames(sigProjMatrix_Pempadj) <- colnames(sigProjMatrix_Pemp)
  rownames(sigProjMatrix_Pempadj) <- rownames(sigProjMatrix_Pemp)

  sigProjMatrix_Pempadj[sigProjMatrix_Pempadj == 0] <- 10^(-300)

  sigProjMatrix_Pempadj <- as.matrix(log10(sigProjMatrix_Pempadj))

  return(list(sigProjMatrix = sigProjMatrix, pVals = sigProjMatrix_Padj, emp_pVals = sigProjMatrix_Pempadj))
}

#' Evaluates the significance of each signature vs. each projection.
#'
#' @importFrom stats p.adjust
#' @param projections Maps projections to their spatial coordinates for each
#' sample
#' @param sigScoresData List of SignatureScores Object, mapping signature
#' names to their value at each coordinate
#' @param metaData data.frame of meta-data for cells
#' @param randomSigData A list with two items:
#'
#'   randomSigs: a list of lists of SignatureScore objects.  Each sub-list
#'     represents permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
#' @return list:
#' \itemize{
#'     \item sigProbMatrix: the matrix of signature-projection consistency scores
#'     \item pVals: pvalues for the scores
#' }
sigsVsProjections <- function(projections, sigScoresData, metaData,
                              randomSigData) {

  sigProjMatrix <- data.frame()
  sigProjMatrix_P <- data.frame()
  sigProjMatrix_Pemp <- data.frame()

  for (proj in projections) {
    weights <- computeKNNWeights(proj, K=round(sqrt(NCOL(proj@pData))))

    svp_n <- sigsVsProjection_n(sigScoresData, randomSigData, weights)
    svp_pcn <- sigsVsProjection_pcn(metaData, weights)
    svp_pcf <- sigsVsProjection_pcf(metaData, weights)

    consistency <- c(svp_n$consistency,
                    svp_pcn$consistency,
                    svp_pcf$consistency)

    consistency <- as.data.frame(consistency)
    colnames(consistency) <- c(proj@name)

    pvals = c(svp_n$pvals,
              svp_pcn$pvals,
              svp_pcf$pvals)

    emp_pvals = c(svp_n$empvals,
                  svp_pcn$pvals,
                  svp_pcf$pvals)

    pvals <- as.data.frame(pvals)
    colnames(pvals) <- c(proj@name)

    emp_pvals <- as.data.frame(emp_pvals)
    colnames(emp_pvals) <- c(proj@name)

    sigProjMatrix <- merge(sigProjMatrix, consistency,
                           by='row.names', all=TRUE)

    rownames(sigProjMatrix) = sigProjMatrix$Row.names
    sigProjMatrix$Row.names <- NULL

    sigProjMatrix_P <- merge(sigProjMatrix_P, pvals,
                             by='row.names', all=TRUE)

    sigProjMatrix_Pemp <- merge(sigProjMatrix_Pemp, emp_pvals,
                                by='row.names', all=TRUE)

    rownames(sigProjMatrix_P) = sigProjMatrix_P$Row.names
    rownames(sigProjMatrix_Pemp) = sigProjMatrix_Pemp$Row.names
    sigProjMatrix_P$Row.names <- NULL
    sigProjMatrix_Pemp$Row.names <- NULL

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

  sigProjMatrix_Pempadj <- matrix(
      p.adjust(
          as.matrix(sigProjMatrix_Pemp),
          method="BH"
      ),
    nrow=nrow(sigProjMatrix_Pemp),
    ncol=ncol(sigProjMatrix_Pemp)
  )
  colnames(sigProjMatrix_Pempadj) <- colnames(sigProjMatrix_Pemp)
  rownames(sigProjMatrix_Pempadj) <- rownames(sigProjMatrix_Pemp)

  sigProjMatrix_Pempadj[sigProjMatrix_Pempadj == 0] <- 10^(-300)

  sigProjMatrix_Pempadj <- as.matrix(log10(sigProjMatrix_Pempadj))

  return(list(sigProjMatrix = sigProjMatrix, pVals = sigProjMatrix_Padj, emp_pVals = sigProjMatrix_Pempadj))
}

#' Evaluates the significance of each numeric signature vs. a
#' single projections weights
#'
#' @importFrom matrixStats colMedians
#' @param sigData a list of SignatureScores objects representing input user
#' signatures
#' @param randomSigData A list with two items:
#'
#'   randomSigs: a list of lists of SignatureScore objects.  Each sub-list
#'     represents permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
#'
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @param cells list of cell names.  Subsets anlysis to provided cells.
#' If omitted, use all cells.
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#'     \item emppvals: empirical pvalues
#' }
sigsVsProjection_n <- function(sigData, randomSigData,
                               weights, cells = NULL){

    # Build a matrix of all non-meta signatures
    if (is.null(cells)){
        extractedScores <- lapply(sigData,
                                  function(sig) sig@scores)
    } else {
        extractedScores <- lapply(sigData,
                                  function(sig) sig@scores[cells])
    }
    sigScoreMatrix <- do.call(cbind, extractedScores)
    sampleLabels <- rownames(sigScoreMatrix)
    sigScoreMatrix <- matrixStats::colRanks(
               sigScoreMatrix, preserveShape = TRUE, ties.method = "average"
    )
    colnames(sigScoreMatrix) <- names(sigData)
    rownames(sigScoreMatrix) <- sampleLabels


    if (!is.null(cells)) {
        weights <- weights[cells, cells]
    }

    N_SAMPLES <- nrow(weights)

    randomSigScores <- randomSigData$randomSigs
    sigAssignments <- randomSigData$sigAssignments

    groupedResults <- lapply(names(randomSigScores), function(group) {

        # Build a matrix of random background signatures for this group
        randomSigGroup <- randomSigScores[[group]]

        if (is.null(cells)){
            extractedScores <- lapply(randomSigGroup,
                                      function(rsig) rsig@scores)
        } else {
            extractedScores <- lapply(randomSigGroup,
                                      function(rsig) rsig@scores[cells])
        }

        randomSigScoreMatrix <- do.call(cbind, extractedScores)

        sampleLabels <- rownames(randomSigScoreMatrix)
        randomSigScoreMatrix <- matrixStats::colRanks(
              randomSigScoreMatrix, preserveShape = TRUE,
              ties.method = "average"
        )

        colnames(randomSigScoreMatrix) <- vapply(
            randomSigGroup, function(x) x@name, ""
        )

        rownames(randomSigScoreMatrix) <- sampleLabels

        groupSigNames <- names(sigAssignments)[sigAssignments == group]

        sigScoreMatrixGroup <- sigScoreMatrix[, groupSigNames]

        # Calculate scores for actual signatures
        neighborhoodPrediction <- as.matrix(weights %*% sigScoreMatrixGroup)
        dissimilarity <- abs(sigScoreMatrixGroup -
                             neighborhoodPrediction)
        medDissimilarity <- colMedians(dissimilarity)
        names(medDissimilarity) <- colnames(dissimilarity)

        # Calculate scores for random signatures
        randomNeighborhoodPrediction <- as.matrix(weights %*% randomSigScoreMatrix)
        randomDissimilarity <- abs(randomSigScoreMatrix -
                                   randomNeighborhoodPrediction)
        randomMedDissimilarity <- colMedians(randomDissimilarity)
        names(randomMedDissimilarity) <- colnames(randomDissimilarity)

        mu <- mean(randomMedDissimilarity)
        sigma <- sd(randomMedDissimilarity)

        #Create CDF function for medDissmilarityPrime and apply CDF function to
        pvals <- pnorm( (medDissimilarity - mu) / sigma)
        consistency <- 1 - (medDissimilarity / N_SAMPLES)

        orderedBg <- sort(randomMedDissimilarity)
        empvals <- vapply(medDissimilarity, function(x) {
                              N <- length(orderedBg)
                              comp <- which(orderedBg < x)
                              if (length(comp) == 0) {
                                  p <- 1 / (N + 1)
                              } else {
                                  p <- (max(comp) + 1) / (N + 1)
                              }
                              return(p)
       }, FUN.VALUE = 0.0)


        return(list(consistency = consistency, pvals = pvals,
                    empvals = empvals))
    })


    consistency <- do.call(c, lapply(groupedResults, function(x) x$consistency))
    pvals <- do.call(c, lapply(groupedResults, function(x) x$pvals))
    empvals <- do.call(c, lapply(groupedResults, function(x) x$empvals))


    return(list(consistency = consistency, pvals = pvals,
                empvals = empvals))
}

#' Evaluates the significance of each meta data numeric signature vs. a
#' single projections weights
#' @importFrom stats pnorm
#' @param metaData data.frame of meta-data for cells
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @param cells list of cell names.  Subsets anlysis to provided cells.
#' If omitted, use all cells.
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#' }
sigsVsProjection_pcn <- function(metaData, weights, cells = NULL){
  # Calculate significance for meta data numerical signatures
  # This is done separately because there are likely to be many repeats
  # (e.g. for a time coordinate)

  if (!is.null(cells)) {
      weights <- weights[cells, cells]
      metaData <- metaData[cells, ]
  }

  N_SAMPLES <- nrow(weights)

  consistency <- numeric()
  pvals <- numeric()

  for (metaName in names(metaData)) {

    scores <- metaData[[metaName]]

    if (is.factor(scores)) {
      next
    }

    sigScores <- rank(scores, ties.method = "average")
    if (all(sigScores == sigScores[1])) {
      consistency[metaName] <- 0
      pvals[metaName] <- 1.0
      next
    }

    sigPredictions <- as.matrix(weights %*% sigScores)
    dissimilarity <- abs(sigScores - sigPredictions)
    medDissimilarity <- median(dissimilarity)

    #Compute a background for numerical signatures
    NUM_REPLICATES <- 3000
    bgValues <- replicate(NUM_REPLICATES, sample(sigScores))
    randomPredictions <- as.matrix(weights %*% bgValues)

    rDissimilarity <- abs(bgValues - randomPredictions)
    randomScores <- as.matrix(apply(rDissimilarity, 2, median))


    mu <- mean(randomScores)
    sigma <- sd(randomScores)

    if (sigma != 0) {
      p_value <- pnorm( (medDissimilarity - mu) / sigma)
    } else {
      p_value <- 1.0
    }

    consistency[metaName] <- 1 - (medDissimilarity / N_SAMPLES)
    pvals[metaName] <- p_value

  }

  return(list(consistency = consistency, pvals = pvals))
}

#' Evaluates the significance of each meta data factor signature vs. a
#' single projections weights
#' @importFrom stats kruskal.test
#' @param metaData data.frame of meta-data for cells
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @param cells list of cell names.  Subsets anlysis to provided cells.
#' If omitted, use all cells.
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#' }
sigsVsProjection_pcf <- function(metaData, weights, cells = NULL){
  # Calculate significance for meta data factor signatures
  # This is done separately because there are likely to be many repeats
  # (e.g. for a time coordinate)

  if (!is.null(cells)) {
      weights <- weights[cells, cells]
      metaData <- metaData[cells, ]
  }

  N_SAMPLES <- nrow(weights)

  consistency <- numeric()
  pvals <- numeric()

  for (sigName in names(metaData)) {

    scores <- metaData[[sigName]]

    if (!is.factor(scores)) {
      next
    }

    ### build one hot matrix for factors

    fValues <- scores
    fLevels <- levels(fValues)
    factorFreq <- matrix(0L, ncol = length(fLevels))
    factorMatrix <- matrix(0L, nrow = N_SAMPLES, ncol = length(fLevels))

    factList <- lapply(fLevels, function(fval) {
      factorMatrixRow <- matrix(0L, nrow = N_SAMPLES, ncol = 1)
      equal_ii <- fValues == fval
      factorMatrixRow[equal_ii] <- 1
      return(list(sum(equal_ii) / length(fValues), factorMatrixRow))
    })
    factorFreq <- lapply(factList, function(x) return(x[[1]]))
    factorMatrix <- matrix(unlist(lapply(factList, function(x) return(x[[2]]))),
                           nrow = N_SAMPLES, ncol = length(fLevels))

    if (1 %in% factorFreq) {
      consistency[sigName] <- 0.0
      pvals[sigName] <- 1.0
      next
    }

    factorPredictions <- as.matrix(weights %*% factorMatrix)

    labels <- apply(factorMatrix, 1, which.max)

    # Get predicted index and corresponding probability
    preds_ii <- apply(factorPredictions, 1, which.max)

    # Calculate the p value for meta data signature
    krList <- list()
    for (k in 1:length(fLevels)) {
        group <- preds_ii[labels == k]
        if (length(group) > 0){
            krList <- c(krList, list(group))
        }
    }

    if (length(krList) > 0){
        krTest <- kruskal.test(krList)

        if (is.na(krTest$p.value)){
            consistency[sigName] <- 0
            pvals[sigName] <- 1.0
        } else {
            consistency[sigName] <- krTest$statistic
            pvals[sigName] <- krTest$p.value
        }
    } else {
        consistency[sigName] <- 0
        pvals[sigName] <- 1.0
    }

  }

  return(list(consistency = consistency, pvals = pvals))
}

#' Clusters signatures according to the rank sum
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom rsvd rsvd
#' @param sigMatrix data.frame of signatures scores, NUM_SAMPLES x NUM_SIGNATURES
#' @param metaData data.frame of meta-data for cells
#' @param pvals the corresponding P-values for each score,
#' NUM_SIGNATURES x NUM_SAMPLES
#' @return a list:
#' \itemize{
#'     \item Computed: a list of clusters of computed signatures
#'     \item Meta: a list of clusters of meta data signatures
#' }
clusterSignatures <- function(sigMatrix, metaData, pvals, clusterMeta) {

  significant <- apply(pvals, 1, function(x) min(x) < -1.3)

  comp_names <- colnames(sigMatrix)
  signif_computed <- significant[comp_names]

  keep_significant <- names(signif_computed)[signif_computed]

  # Cluster computed signatures and precomputed signatures separately
  computedSigMatrix <- sigMatrix[, keep_significant, drop = FALSE]

  compcls <- list()
  maxcls <- 1
  if (ncol(computedSigMatrix) > 1) {

    if (nrow(computedSigMatrix) > 5000) {
        cSM_sub <- computedSigMatrix[sample(nrow(computedSigMatrix), 5000), ]
    } else {
        cSM_sub <- computedSigMatrix
    }

    r <- colRanks(
            as.matrix(cSM_sub),
            ties.method = "average",
            preserveShape = TRUE
         )

    r <- t(r)

    mbic <- mclustBIC(r, modelNames = "VVI")
    compkm <- Mclust(r, x = mbic)

    compcls <- as.list(compkm$classification)
    names(compcls) <- colnames(computedSigMatrix)
    compcls <- compcls[order(unlist(compcls), decreasing = FALSE)]

    maxcls <- max(unlist(compcls))
  }

  compcls[names(signif_computed)[!signif_computed]] <- maxcls + 1

  # Cluster meta data signatures by grouping % signatures with their
  # factors if they exist
  metacls <- list()
  if (!is.null(metaData) && ncol(metaData) > 0) {
      metaSigsToCluster <- colnames(metaData)
      metacls[metaSigsToCluster] <- 1

      if (clusterMeta) {
          pc_i <- 2
          for (name in metaSigsToCluster) {
              if (is.factor(metaData[[name]])) {
                  metacls[[name]] <- pc_i
                  sigLevels <- levels(metaData[[name]])
                  for (level in sigLevels) {
                      key <- paste(name, level, sep = "_")
                      if (key %in% names(metacls)) {
                          metacls[[key]] <- pc_i
                      }
                  }
                  pc_i <- pc_i + 1
              }
          }
      }
  }

  return(list(Computed = compcls, Meta = metacls))
}
