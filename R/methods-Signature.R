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
#'
#' @param sigDict Named list of signs for each gene in the signature
#' @param name Name of the signature
#' @param source File from which this signature was read from
#' @param metaData Metadata pertinent to signature
#' @return Signature object
Signature <- function(sigDict, name, source, metaData="") {
  if (missing(sigDict)) {
    stop("Missing sigDict information.")
  } else if (missing(name)) {
    stop("Missing signature name.")
  } else if (missing(source)) {
    stop("Missing source file.")
  }

  names(sigDict) <- toupper(names(sigDict))

  .Object <- new("Signature", sigDict = sigDict, name = name,
                 source = source, metaData = metaData)

  return(.Object)
}

#' Create a user-defined gene signature
#'
#' Typical usage of VISION involes providing a location to a ".gmt" signature
#' file from which VISION automatically creates Signature objects.  However,
#' using the createGeneSignature method, users may programmatically define
#' signatures from within R.
#'
#' @param name the name of the signature
#' @param sigData a named vector where the names correspond to genes in the
#' data and the values are either `1.0` for up-regulated (or positive) genes,
#' and `-1.0` for down-regulated (negtive) genes. For an unsigned signature
#' use 1.0 for all values.
#' @param metadata metadata that is relevent to the signature. [Default:NULL]
#' @export
#' @return a Signature object
#' @examples
#' \dontrun{
#' sig1 <- createGeneSignature(
#'            name = "CD8 Markers",
#'            sigData = c(CD8A=1, CD8B=1, GZMK=1, GZMB=1,
#'                        GZMH=1, GZMA=1, GNLY=1, DUSP2=1,
#'                        EOMES=1, TBX21=1, PRMD1=1, PRF1=1,
#'                        IFNG=1)
#'        )
#'
#' cc_sigs <- "path/to/cell_cycle.gmt"
#'
#' sigs <- c(sig1, cc_sigs)
#'
#' vis <- Vision(data = expMat, signatures = sigs)
#' }
createGeneSignature <- function(name, sigData, metadata="") {
    return(
        Signature(sigData, name, source = "user-defined", metaData = metadata)
    )
}

#' Generate random signatures for a null distribution by permuting the data
#'
#' @importFrom stats runif
#' @param eData the data to use for the permutations
#' @param sigData list of signature objects
#' random signature sizes
#' @param num the number of signatures to generate
#' @return A list with two items:
#'
#'   randomSigs: a list of lists of Signature objects.  Each sub-list represents
#'     permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
generatePermutationNull <- function(eData, sigData, num) {

  exp_genes <- rownames(eData)

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
      clusters <- as.factor(seq_len(nrow(sigVars)))
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
  rownames(centers) <- as.character(seq_len(n_components))

  # undo the log scaling
  centers[, "sigSize"] <- round(10 ** centers[, "sigSize"])

  message("Creating ", nrow(centers),
          " background signature groups with the following parameters:")
  print(centers) # How do I do this with 'message'??
  message("  signatures per group: ", num)

  randomSigs <- list()
  randomSigAssignments <- character()

  for (cluster_i in rownames(centers)) {

    size <- centers[cluster_i, "sigSize"]
    balance <- centers[cluster_i, "sigBalance"]

    for (j in 1:num) {
      newSigGenes <- sample(exp_genes, min(size, length(exp_genes)))

      upGenes <- floor(balance * size)
      remainder <- (balance * size) %% 1
      if (runif(1, 0, 1) < remainder){
          upGenes <- upGenes + 1
      }
      newSigSigns <- c(rep(1, upGenes), rep(-1, size - upGenes))

      names(newSigSigns) <- newSigGenes
      newSig <- Signature(
          newSigSigns, paste0("RANDOM_BG_", cluster_i, "_", j), "x")
      randomSigs[[newSig@name]] <- newSig
      randomSigAssignments[[newSig@name]] <- cluster_i
    }

  }

  randomSigAssignments <- as.factor(randomSigAssignments)

  return(
      list(
          randomSigs = randomSigs,
          sigAssignments = clusters,
          randomSigAssignments = randomSigAssignments
          )
  )
}

#' Evaluates the significance of each signature in each cluster
#'
#' @importFrom stats p.adjust
#' @param weights output of computeKNNWeights
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
#' @param normExpr NormData object used to generate random background sigs
#' @return list:
#' \itemize{
#'     \item sigProbMatrix: the vector of consistency z-scores
#'     \item pVals: pvalues for the scores
#'     \item emp_pVals: pvalues for the scores
#' }
sigConsistencyScores <- function(weights, sigScoresData,
                                 metaData, randomSigData,
                                 normExpr) {

    signatureNames <- c(colnames(sigScoresData), colnames(metaData))

    svp_n <- sigsVsProjection_n(
        sigScoresData, randomSigData, normExpr, weights)
    svp_pcn <- sigsVsProjection_pcn(metaData, weights)
    svp_pcf <- sigsVsProjection_pcf(metaData, weights)

    consistency <- c(svp_n$consistency,
                    svp_pcn$consistency,
                    svp_pcf$consistency)

    pvals <- c(svp_n$pvals,
              svp_pcn$pvals,
              svp_pcf$pvals)

    consistency <- consistency[signatureNames]
    pvals <- pvals[signatureNames]

    # FDR-correct and log-transform p-values
    if (length(pvals) > 1){
        fdr <- p.adjust(pvals, method = "BH")
    } else {
        fdr <- numeric(0)
    }

    df <- data.frame(
        C = consistency,
        pValue = pvals,
        FDR = fdr
    )

    return(df)
}


#' Evaluates the significance of each protein
#'
#' @importFrom stats p.adjust
#' @param weights output of computeKNNWeights
#' @param proteinData numeric matrix of protein abundance
#' size is cells x proteins
#' @return dataframe with columns "C", "pValue", and "FDR"
fbConsistencyScores <- function(weights, proteinData) {

    proteinNames <- colnames(proteinData)
    svp_pcn <- sigsVsProjection_pcn(proteinData, weights, computePval = FALSE)

    consistency <- svp_pcn$consistency

    pvals <- svp_pcn$pvals

    consistency <- consistency[proteinNames]
    pvals <- pvals[proteinNames]

    result <- data.frame(
        C = consistency, pValue = pvals
    )
    result$FDR <- p.adjust(result$pValue, method = "BH")

    return(result)
}


#' Evaluates the significance of each numeric signature vs. a
#' single projections weights
#'
#' @importFrom matrixStats colMedians
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats setNames
#' @param sigScores numeric matrix of signature scores
#' size is cells x signatures
#' @param randomSigData A list with three items:
#'
#'   randomSigs: a list of signature score matrices.  Each list item
#'     represents permutation signatures generated for a specific size/balance,
#'     and is a numeric matrix of size cells X signatures
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
#'
#'   randomSigAssignments: named factor vector assigning signatures to random background
#'     groups
#'
#' obtained from calling generatePermutationNull
#' @param normExpr NormData object used to generate random background sigs
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @param cells list of cell names.  Subsets anlysis to provided cells.
#' If omitted, use all cells.
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#'     \item emppvals: empirical pvalues
#' }
sigsVsProjection_n <- function(sigScores, randomSigData,
                               normExpr, weights, cells = NULL){

    if (length(sigScores) == 0){
        return(list(consistency = numeric(), pvals = numeric(),
                empvals = numeric()))
    }

    # Build a matrix of all non-meta signatures
    sigScoreMatrix <- sigScores

    if (!is.null(cells)){
        sigScoreMatrix <- sigScoreMatrix[cells, , drop = FALSE]
    }

    if (!is.null(cells)) {
        weights$indices <- weights$indices[cells, , drop = FALSE]
        weights$weights <- weights$weights[cells, , drop = FALSE]
    }

    gearyFG <- pbmclapply(seq_len(ncol(sigScoreMatrix)), function(ii) {
        sigScoreMatrixGroup <- sigScoreMatrix[, ii, drop = FALSE]
        sigScoreMatrixGroup <- matrixStats::colRanks(
                   sigScoreMatrixGroup, preserveShape = TRUE, ties.method = "average"
        )
        geary_c <- geary_sig_v_proj(sigScoreMatrixGroup, weights$indices, weights$weights)
        return(geary_c)
    }, mc.preschedule = FALSE)
    gearyFG <- unlist(gearyFG)
    names(gearyFG) <- colnames(sigScoreMatrix)

    randomSigs <- randomSigData$randomSigs
    sigAssignments <- randomSigData$sigAssignments
    randomSigAssignments <- randomSigData$randomSigAssignments

    randomSigBatches <- batchify(randomSigs, 10)

    gearyBG <- pbmclapply(randomSigBatches, function(randomSigSubset){
        randomSigScores <- t(innerEvalSignatureBatchNorm(normExpr, randomSigSubset))
        randomSigScores <- matrixStats::colRanks(
            randomSigScores, preserveShape = TRUE, ties.method = "average"
        )
        geary_c <- geary_sig_v_proj(randomSigScores, weights$indices, weights$weights)
        return(geary_c)
    }, mc.preschedule = FALSE)
    gearyBG <- setNames(unlist(gearyBG), names(randomSigs))

    # Compute fg with bg
    pvals <- lapply(levels(sigAssignments), function(level){
        fgNames <- names(sigAssignments)[sigAssignments == level]
        fg <- gearyFG[fgNames]

        bgNames <- names(randomSigAssignments)[randomSigAssignments == level]
        bg <- gearyBG[bgNames]

        N <- length(bg)
        empvals <- vapply(fg, function(x) {
            comp <- sum(bg < x)
            p <- (comp + 1) / (N + 1)
            return(p)
        }, FUN.VALUE = 0.0)
    })
    pvals <- unlist(pvals)

    return(list(consistency = 1 - gearyFG, pvals = pvals))
}

#' Evaluates the significance of each meta data numeric signature vs. a
#' single projections weights
#' @importFrom stats pnorm
#' @importFrom stats setNames
#' @importFrom pbmcapply pbmclapply
#' @param metaData data.frame of meta-data for cells
#' @param weights numeric matrix of dimension N_SAMPLES x N_SAMPLES
#' @param cells list of cell names.  Subsets anlysis to provided cells.
#' If omitted, use all cells.
#' @param computePval bool whether to compute p-values. Set to FALSE to
#' save computation time for large numbers of numeric signatures. If set
#' to FALSE, all p-values are reported as 1.0
#' @return list:
#' \itemize{
#'     \item consistency: consistency scores
#'     \item pvals: pvalues
#' }
sigsVsProjection_pcn <- function(metaData, weights, cells = NULL, computePval = TRUE){
    # Calculate significance for meta data numerical signatures
    # This is done separately because there are likely to be many repeats
    # (e.g. for a time coordinate)

    if (!is.null(cells)) {
        weights$indices <- weights$indices[cells, , drop = FALSE]
        weights$weights <- weights$weights[cells, , drop = FALSE]
        metaData <- metaData[cells, , drop = FALSE]
    }

    if (is.data.frame(metaData)){
        numericMeta <- vapply(names(metaData),
                               function(metaName) {
                                   is.numeric(metaData[, metaName])
                               }, FUN.VALUE = TRUE)

        numericMetaNames <- names(numericMeta)[numericMeta]
    } else { # If it's a matrix, assume all numeric
        numericMetaNames <- colnames(metaData)
    }

    if (length(numericMetaNames) == 0) {
        return(list(consistency = c(), pvals = c()))
    }

    numericMetaRanks <- matrixStats::colRanks(
        as.matrix(metaData[, numericMetaNames, drop = FALSE]),
        preserveShape = TRUE, ties.method = "average"
    )
    colnames(numericMetaRanks) <- numericMetaNames

    NUM_REPLICATES <- 3000

    # Generate list of jobs
    type <- character()
    col <- numeric()
    for (i in colnames(numericMetaRanks)) {
        type <- c(type, "fg")
        col <- c(col, i)
        if (computePval){
            type <- c(type, replicate(NUM_REPLICATES, "bg"))
            col <- c(col, replicate(NUM_REPLICATES, i))
        }
    }
    jobs <- data.frame(type = type, col = col, stringsAsFactors = FALSE)

    results <- pbmclapply(seq_len(nrow(jobs)), function(i) {

        type <- jobs[i, "type"]
        col <- jobs[i, "col"]


        sigScores <- numericMetaRanks[, col, drop = FALSE]

        if (any(is.na(sigScores))){
            return(0.0)
        }

        if (all(sigScores == sigScores[1])) {
            return(0.0)
        }

        if (type == "bg"){
            sigScores <- matrix(sample(sigScores), ncol = 1)
            colnames(sigScores) <- colnames(numericMetaRanks)[col]
        }

        geary_c <- geary_sig_v_proj(
            sigScores, weights$indices, weights$weights)

        return(1 - geary_c)

    }, mc.preschedule = FALSE)

    jobs[["C"]] <- unlist(results)

    results <- jobs[jobs$type == "fg", c("col", "C")]
    rownames(results) <- results$col
    results$pval <- 1.0

    if (computePval) {
        for (varname in rownames(results)) {
            bg <- jobs[(jobs$col == varname) & (jobs$type == "bg"), ]
            comp <- sum(bg$C >= results[varname, "C"])
            pval <- (comp + 1) / (NUM_REPLICATES + 1)
            results[varname, "pval"] <- pval
        }
    }

    consistency <- setNames(results$C, rownames(results))
    pvals <- setNames(results$pval, rownames(results))

    return(list(consistency = consistency, pvals = pvals))
}

#' Evaluates the significance of each meta data factor signature vs. a
#' single projections weights
#' @importFrom stats chisq.test
#' @importFrom pbmcapply pbmclapply
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

  results <- pbmclapply(factorMeta, function(metaName) {

    scores <- metaData[[metaName]]

    ### build one hot matrix for factors

    fValues <- droplevels(scores)
    fLevels <- levels(fValues)

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

    contingency_rows <- lapply(seq(fLevels), function(k){
        group <- factorPredictions[labels == k, , drop = FALSE]
        return(colSums(group))
    })

    contingency <- do.call(rbind, contingency_rows)

    if (nrow(contingency) > 1){

        # Drop cols where all 0
        # Should only happen in very rare circumstances
        # (since neighborhoods aren't symmetric)
        good_cols <- colSums(contingency) > 0
        contingency <- contingency[, good_cols]

        suppressWarnings({
            chsqResults <- chisq.test(contingency)
        })

        if (is.na(chsqResults$p.value)){
            c_score <- 0.0
            pval <- 1.0
        } else {
            # for the c_score, compute the 1-V (cramers V)
            n <- sum(contingency)
            V <- sqrt(chsqResults$statistic / n /
                      min(nrow(contingency) - 1, ncol(contingency) - 1)
                  )
            c_score <- V
            pval <- chsqResults$p.value
        }
    } else {
        c_score <- 0.0
        pval <- 1.0
    }

    return(list(consistency = c_score, pval = pval))

  })

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
#'
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom rsvd rsvd
#' @param sigMatrix matrix of signatures scores, NUM_SAMPLES x NUM_SIGNATURES
#' @param metaData data.frame of meta-data for cells
#' @param autocorrelation data.frame autocorrelation results
#' @param clusterMeta bool whether or not to cluster meta-data (used when pool=TRUE)
#' @return a list:
#' \itemize{
#'     \item Computed: a list of clusters of computed signatures
#'     \item Meta: a list of clusters of meta data signatures
#' }
clusterSignatures <- function(sigMatrix, metaData, autocorrelation, clusterMeta) {

  significant <- autocorrelation$FDR < .05

  # Additionally, threshold on the Geary's C' itself
  large <- autocorrelation$C > 0.2
  significant <- significant & large
  names(significant) <- rownames(autocorrelation)

  meta_n <- vapply(names(metaData), function(metaName) {
			is.numeric(metaData[, metaName]) && !any(is.na(metaData[, metaName]))
		}, FUN.VALUE = TRUE)

  meta_n <- metaData[, meta_n, drop = FALSE]
  sigMatrix <- cbind(sigMatrix, meta_n)

  comp_names <- colnames(sigMatrix)
  signif_computed <- significant[comp_names]

  keep_significant <- names(signif_computed)[signif_computed]

  # Cluster computed signatures and precomputed signatures separately
  computedSigMatrix <- sigMatrix[, keep_significant, drop = FALSE]

  compcls <- list()
  maxcls <- 1
  if (ncol(computedSigMatrix) > 5) {

    # z-normalize each signature vector before clustering
    computedSigMatrix <- colNormalization(as.matrix(computedSigMatrix))


    if (nrow(computedSigMatrix) > 5000) {
        cSM_sub <- computedSigMatrix[sample(nrow(computedSigMatrix), 5000), ]
    } else {
        cSM_sub <- computedSigMatrix
    }

    r <- t(as.matrix(cSM_sub))

    mbic <- mclustBIC(r, G=1:15, modelNames = "EII")
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
