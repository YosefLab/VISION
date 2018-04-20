#' Processes signatures on input
#'
#' Removes genes that don't match a gene in the expression matrix
#' Drops signatures with less than `minSignatureGenes` matching genes
#'
#' @param sigData list of signature objects
#' @param expressionGenes list of gene identifiers in the expression matrix
#' @param minSignatureGenes minimum number of genes a signature must match
#' in the expression matrix in order to be retained
#' @return processedSigData list of signature objects
processSignatures <- function(sigData, expressionGenes, minSignatureGenes){
    out <- lapply(sigData, function(sig){
        validGenes <- names(sig@sigDict) %in% expressionGenes
        sig@sigDict <- sig@sigDict[validGenes]
        return(sig)
    })

    names(out) <- names(sigData)

    validSigs <- vapply(out,
                        function(sig) length(sig@sigDict) >= minSignatureGenes,
                        FUN.VALUE = FALSE)

    out <- out[validSigs]
    return(out)
}

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
                            normExpr, object@weights)

    # randomSigScores is matrix of cells x signatures
    # need to make a list of matrices, one for each signature group
    # by background group number
    randomSigScoresGroups <- list()
    for (groupName in names(randomSigs)) {
        sigsInGroup <- names(randomSigs[[groupName]])
        randomSigScoresGroups[[groupName]] <- randomSigScores[, sigsInGroup, drop = FALSE]
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
                        return(length(s@sigDict))
                    }
                 , 1)

  sigSize <- log10(sigSize)
  sigBalance <- vapply(sigData, function(s) {
                           positive <- sum(s@sigDict >= 0)
                           balance <- positive / length(s@sigDict)
                           return(balance)
    }, 1)

  sigBalance[sigBalance < 0.5] <- 1 - sigBalance[sigBalance < 0.5]

  sigVars <- cbind(sigSize, sigBalance)

  n_components <- 5 # TODO: choose number of components better

  if (nrow(sigVars) <= n_components){
      n_components <- nrow(sigVars)
      centers <- sigVars
      clusters <- as.factor(seq(1:nrow(sigVars)))
      names(clusters) <- rownames(sigVars)
  } else {

      if (nrow(unique(sigVars)) <= n_components){
          n_components <- nrow(unique(sigVars))
      }

      km <- kmeans(sigVars, n_components)
      centers <- km$centers

      levels <- as.character(seq(n_components))
      clusters <- factor(km$cluster, levels = levels)
  }

  # Re-order the centers
  row_i <- order(centers[, "sigSize"], centers[, "sigBalance"])

  centers <- centers[row_i, , drop = FALSE]
  levels(clusters) <- as.character(order(row_i))
  rownames(centers) <- as.character(seq(n_components))

  # undo the log scaling
  centers[, "sigSize"] <- round(10 ** centers[, "sigSize"])

  message("Creating ", nrow(centers),
          " background signature groups with the following parameters:")
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

#' Evaluates the significance of each signature in each cluster
#'
#' @importFrom stats p.adjust
#' @param latentSpace numeric matrix N_Cells x N_Components
#' @param sigScoresData numeric matrix of signature scores
#' size is cells x signatures
#' @param metaData data.frame of meta-data for cells
#' @param randomSigData A list with two items:
#'
#'   randomSigs: a list of signature score matrices.  Each list item
#'     represents permutation signatures generated for a specific size/balance,
#'     and is a numeric matrix of size cells X signatures
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
#' @param fdrCorrect logical, whether or not to FDR correct p-values
#' @return list:
#' \itemize{
#'     \item sigProbMatrix: the vector of consistency z-scores
#'     \item pVals: pvalues for the scores
#'     \item emp_pVals: pvalues for the scores
#' }
sigConsistencyScores <- function(latentSpace, sigScoresData,
                                      metaData, randomSigData,
                                      fdrCorrect = TRUE) {

    signatureNames <- c(colnames(sigScoresData), colnames(metaData))

    weights <- computeKNNWeights(latentSpace)

    svp_n <- sigsVsProjection_n(sigScoresData,
                                randomSigData, weights)
    svp_pcn <- sigsVsProjection_pcn(metaData, weights)
    svp_pcf <- sigsVsProjection_pcf(metaData, weights)

    consistency <- c(svp_n$consistency,
                    svp_pcn$consistency,
                    svp_pcf$consistency)

    pvals <- c(svp_n$pvals,
              svp_pcn$pvals,
              svp_pcf$pvals)

    emp_pvals <- c(svp_n$empvals,
                  svp_pcn$pvals,
                  svp_pcf$pvals)

    consistency <- consistency[signatureNames]
    pvals <- pvals[signatureNames]
    emp_pvals <- emp_pvals[signatureNames]

    consistency <- as.matrix(consistency)
    pvals <- as.matrix(pvals)
    emp_pvals <- as.matrix(emp_pvals)

    colnames(consistency) <- c("Consistency")
    colnames(pvals) <- c("Consistency")
    colnames(emp_pvals) <- c("Consistency")

    # Cast to 1-column matrices so its consistent with other ProjectionData outputs

    # FDR-correct and log-transform p-values
    if (fdrCorrect){
        pvals <- apply(pvals, MARGIN = 2,
                                    FUN = p.adjust, method = "BH")
        pvals[pvals == 0] <- 10 ^ (-300)
        pvals <- log10(pvals)

        emp_pvals <- apply(emp_pvals, MARGIN = 2,
                                       FUN = p.adjust, method = "BH")
        emp_pvals[emp_pvals == 0] <- 10 ^ (-300)
        emp_pvals <- log10(emp_pvals)
    }

    return(list(sigProjMatrix = consistency,
                pVals = pvals,
                emp_pVals = emp_pvals)
    )
}


#' Evaluates the significance of each numeric signature vs. a
#' single projections weights
#'
#' @importFrom matrixStats colMedians
#' @importFrom parallel mclapply
#' @param sigData numeric matrix of signature scores
#' size is cells x signatures
#' @param randomSigData A list with two items:
#'
#'   randomSigs: a list of signature score matrices.  Each list item
#'     represents permutation signatures generated for a specific size/balance,
#'     and is a numeric matrix of size cells X signatures
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
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
    sigScoreMatrix <- sigData
    if (!is.null(cells)){
        sigScoreMatrix <- sigScoreMatrix[cells, , drop = FALSE]
    }

    rowLabels <- rownames(sigScoreMatrix)
    colLabels <- colnames(sigScoreMatrix)
    sigScoreMatrix <- matrixStats::colRanks(
               sigScoreMatrix, preserveShape = TRUE, ties.method = "average"
    )
    colnames(sigScoreMatrix) <- colLabels
    rownames(sigScoreMatrix) <- rowLabels


    if (!is.null(cells)) {
        weights$indices <- weights$indices[cells, , drop = FALSE]
        weights$weights <- weights$weights[cells, , drop = FALSE]
    }

    N_SAMPLES <- nrow(weights$indices)

    randomSigScores <- randomSigData$randomSigs
    sigAssignments <- randomSigData$sigAssignments

    availableCores <- max(parallel::detectCores() - 1, 1)
    groupedResults <- mclapply(names(randomSigScores), function(group) {

        # Build a matrix of random background signatures for this group
        randomSigScoreMatrix <- randomSigScores[[group]]

        if (!is.null(cells)){
            randomSigScoreMatrix <- randomSigScoreMatrix[cells, , drop = FALSE]
        }

        rowLabels <- rownames(randomSigScoreMatrix)
        colLabels <- colnames(randomSigScoreMatrix)
        randomSigScoreMatrix <- matrixStats::colRanks(
              randomSigScoreMatrix, preserveShape = TRUE,
              ties.method = "average"
        )

        colnames(randomSigScoreMatrix) <- colLabels
        rownames(randomSigScoreMatrix) <- rowLabels

        groupSigNames <- names(sigAssignments)[sigAssignments == group]
        groupSigNames <- intersect(groupSigNames, colnames(sigScoreMatrix))

        sigScoreMatrixGroup <- sigScoreMatrix[, groupSigNames, drop = FALSE]

        # Calculate scores for actual signatures
        medDissimilarity <- geary_sig_v_proj(sigScoreMatrixGroup, weights$indices, weights$weights)

        # Calculate scores for random signatures
        randomMedDissimilarity <- geary_sig_v_proj(randomSigScoreMatrix, weights$indices, weights$weights)

        mu <- mean(randomMedDissimilarity)
        sigma <- sd(randomMedDissimilarity)

        #Create CDF function for medDissmilarityPrime and apply CDF function to
        pvals <- pnorm( (medDissimilarity - mu) / sigma)
        consistency <- (medDissimilarity - mu) / sigma * -1 # Make z-score, higher better

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

    }, mc.cores = max(min(availableCores, length(randomSigScores)), 1))


    consistency <- do.call(c, lapply(groupedResults, function(x) x$consistency))
    pvals <- do.call(c, lapply(groupedResults, function(x) x$pvals))
    empvals <- do.call(c, lapply(groupedResults, function(x) x$empvals))


    return(list(consistency = consistency, pvals = pvals,
                empvals = empvals))
}

#' Evaluates the significance of each meta data numeric signature vs. a
#' single projections weights
#' @importFrom stats pnorm
#' @importFrom parallel mclapply
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
      weights$indices <- weights$indices[cells, , drop = FALSE]
      weights$weights <- weights$weights[cells, , drop = FALSE]
      metaData <- metaData[cells, , drop = FALSE]
  }

  N_SAMPLES <- nrow(weights$indices)

  numericMeta <- vapply(names(metaData),
                         function(metaName) {
                             is.numeric(metaData[, metaName])
                         }, FUN.VALUE = TRUE)

  numericMeta <- names(numericMeta)[numericMeta]

  results <- mclapply(numericMeta, function(metaName) {

    scores <- metaData[[metaName]]

    sigScores <- rank(scores, ties.method = "average")
    if (all(sigScores == sigScores[1])) {
      return(list(consistency = 0.0, pval = 1.0))
    }

    sigScores <- matrix(sigScores, ncol=1)
    rownames(sigScores) <- rownames(metaData)
    colnames(sigScores) <- metaName

    medDissimilarity <- geary_sig_v_proj(sigScores,
                                         weights$indices,
                                         weights$weights)

    #Compute a background for numerical signatures
    NUM_REPLICATES <- 3000
    bgValues <- replicate(NUM_REPLICATES, sample(sigScores))
    rownames(bgValues) <- rownames(sigScores)
    randomScores <- geary_sig_v_proj(bgValues,
                                         weights$indices,
                                         weights$weights)

    mu <- mean(randomScores)
    sigma <- sd(randomScores)

    orderedBg <- sort(as.numeric(randomScores))
    N <- length(orderedBg)
    comp <- which(orderedBg < medDissimilarity)

    if (length(comp) == 0) {
        pval <- 1 / (N + 1)
    } else {
        pval <- (max(comp) + 1) / (N + 1)
    }

    c_score <- (medDissimilarity - mu) / sigma * -1 # make z-score, higher better

    return(list(consistency = c_score, pval = pval))

  }, mc.cores = max(min(10, length(numericMeta)), 1))

  names(results) <- numericMeta

  consistency <- vapply(results, function(x) x$consistency,
                        FUN.VALUE = 0.0)

  pvals <- vapply(results, function(x) x$pval,
                        FUN.VALUE = 0.0)

  return(list(consistency = consistency, pvals = pvals))
}

#' Evaluates the significance of each meta data factor signature vs. a
#' single projections weights
#' @importFrom stats kruskal.test
#' @importFrom parallel mclapply
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

  # load weights into a sparse matrix
  tnn <- t(weights$indices)
  j <- as.numeric(tnn)
  i <- as.numeric(col(tnn))
  vals <- as.numeric(t(weights$weights))
  dims <- c(nrow(weights$indices), nrow(weights$indices))
  dimnames <- list(rownames(weights$indices), rownames(weights$indices))

  weights <- sparseMatrix(i = i, j = j, x = vals,
                          dims = dims, dimnames = dimnames)


  if (!is.null(cells)) {
      weights <- weights[cells, cells]
      metaData <- metaData[cells, , drop = FALSE]
  }

  N_SAMPLES <- nrow(weights)

  factorMeta <- vapply(names(metaData),
                         function(metaName) {
                             is.factor(metaData[, metaName])
                         }, FUN.VALUE = TRUE)

  factorMeta <- names(factorMeta)[factorMeta]

  results <- mclapply(factorMeta, function(metaName) {

    scores <- metaData[[metaName]]

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
        return(list(consistency = 0.0, pval = 1.0))
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
            c_score <- 0
            pval <- 1.0
        } else {
            # for the c_score, approximate the z-score using the chi-square dist
            df <- length(krList) - 1
            c_score <- (krTest$statistic - df ) / sqrt(2 * df)
            pval <- krTest$p.value
        }
    } else {
        c_score <- 0
        pval <- 1.0
    }

    return(list(consistency = c_score, pval = pval))

  }, mc.cores = max(min(10, length(factorMeta)), 1))

  names(results) <- factorMeta

  consistency <- vapply(results, function(x) x$consistency,
                        FUN.VALUE = 0.0)

  pvals <- vapply(results, function(x) x$pval,
                        FUN.VALUE = 0.0)

  return(list(consistency = consistency, pvals = pvals))
}


#' Evaluates values vs coordinates using the Geary C
#'
#' @param values numeric matrix of dimension N_SAMPLES x N_SIGNATURES
#' @param indices numeric matrix of dimension N_SAMPLES x N_NEIGHBORS
#' @param weights numeric matrix of dimension N_SAMPLES x N_NEIGHBORS
#' @return gearyC test statistic values for each signature.  Numeric vector
#' of length N_SIGNATURES
geary_sig_v_proj <- function(values, indices, weights){

    if (nrow(values) != nrow(indices)){
        stop("`values` and `indices` must have same row count")
    }

    if (nrow(weights) != nrow(indices)){
        stop("`weights` and `indices` must have same row count")
    }

    if (ncol(weights) != ncol(indices)){
        stop("`weights` and `indices` must have same column count")
    }

    result <- geary_sparse_all(t(values), indices, weights)
    names(result) <- colnames(values)
    return(result)
}

#' Clusters signatures according to the rank sum
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom rsvd rsvd
#' @param sigMatrix matrix of signatures scores, NUM_SAMPLES x NUM_SIGNATURES
#' @param metaData data.frame of meta-data for cells
#' @param pvals the corresponding P-values for each signature
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
  if (ncol(computedSigMatrix) > 5) {

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
  } else {
      compcls <- as.list(rep(maxcls, ncol(computedSigMatrix)))
      names(compcls) <- colnames(computedSigMatrix)
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
