#' Creates clustering of the cells
#'
#' Results of this are stored as a new variabe in the object's metaData
#'
#' @param object the VISION object for which to cluster the cells
#' @return the VISION object modifed as described above
clusterCells <- function(object) {

    message("Clustering cells...", appendLF = FALSE)

    res <- object@LatentSpace

    n_workers <- getOption("mc.cores")
    n_workers <- if (is.null(n_workers)) 2 else n_workers

    K <- min(object@params$numNeighbors, 30)
    kn <- ball_tree_knn(res,
                        K,
                        n_workers)

    cl <- louvainCluster(kn, res)

    names(cl) <- paste('Cluster', seq(length(cl)))

    # cl is list of character vector
    cluster_variable <- "VISION_Clusters"
    metaData <- object@metaData

    metaData[cluster_variable] <- factor(levels = names(cl))

    for (cluster in names(cl)) {
        metaData[cl[[cluster]], cluster_variable] <- cluster
    }

    object@metaData <- metaData

    message("completed\n")

    return(object)

}


#' create micro-clusters that reduce noise and complexity while maintaining
#' the overall signal in the data
#' @param object the VISION object for which to cluster the cells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @return the VISION with pooled cells
poolCells <- function(object, cellsPerPartition = NULL) {

    if (!is.null(cellsPerPartition)){
        object@params$micropooling$cellsPerPartition <- cellsPerPartition
    }

    if (length(object@Pools) == 0) {
        message(paste(
          "Performing micro-pooling on",
          ncol(object@exprData),
          "cells with a target pool size of",
          object@params$micropooling$cellsPerPartition
        ))
    } else {
      message("Performing micro-pooling using precomputed pools")
    }

    preserve_clusters <- NULL

    if (length(object@Pools) == 0) {
        pools <- applyMicroClustering(
            object@exprData,
            cellsPerPartition = object@params$micropooling$cellsPerPartition,
            filterInput = object@params$latentSpace$projectionGenes,
            filterThreshold = object@params$latentSpace$threshold,
            latentSpace = object@LatentSpace,
            K = object@params$numNeighbors)

        object@Pools <- pools
    }

    message("    Aggregating data using assigned pools...", appendLF = FALSE)
    pooled_cells <- poolMatrixCols(object@exprData, object@Pools)
    object@exprData <- pooled_cells

    if (hasUnnormalizedData(object)) {
        pooled_unnorm <- poolMatrixCols(object@unnormalizedData, object@Pools)
        object@unnormalizedData <- pooled_unnorm
    }

    if (!all(dim(object@LatentSpace) == c(1, 1))) {
        pooled_latent <- poolMatrixRows(object@LatentSpace, object@Pools)
        object@LatentSpace <- pooled_latent
    }

    if (hasFeatureBarcodeData(object)) {
        pooled_fbc <- poolMatrixRows(object@featureBarcodeData, object@Pools)
        object@featureBarcodeData <- pooled_fbc
    }

    poolMeta <- poolMetaData(object@metaData, object@Pools)
    object@metaData <- poolMeta

    if (length(object@Projections) > 0){

        newProjections <- lapply(
            object@Projections,
            function(proj) {
                new_coords <- poolMatrixRows(proj, object@Pools)
                return(new_coords)
            })

        names(newProjections) <- names(object@Projections)

        object@Projections <- newProjections

    }

    message("completed\n")

    return(object)
}

#' filter data accourding to the provided filters
#' @param object the VISION object
#' @param threshold threshold to apply for the threshold filter
#' @param projection_genes either a list of genes or a method to select genes
#' @return the VISION object, populated with filtered data
filterData <- function(object,
                       threshold = NULL,
                       projection_genes = NULL) {

    if (!is.null(threshold)){
        object@params$latentSpace$threshold <- threshold
    }

    if (!is.null(projection_genes)){
        object@params$latentSpace$projectionGenes <- projection_genes
    }

    message("Determining projection genes...")

    if (length(object@params$latentSpace$projectionGenes) == 1){

        exprData <- matLog2(object@exprData)
        projection_genes <- applyFilters(
                    exprData,
                    object@params$latentSpace$threshold,
                    object@params$latentSpace$projectionGenes)

        if (length(projection_genes) == 0){
            stop(
                sprintf("Filtering with (projection_genes=\"%s\", threshold=%i) results in 0 genes\n  Set a lower threshold and re-run",
                    object@params$latentSpace$projectionGenes, object@params$latentSpace$threshold)
                )
        }
    } else {
        projection_genes <- intersect(
            object@params$latentSpace$projectionGenes, rownames(object@exprData))

        if (length(projection_genes) == 0){
            stop("Supplied list of genes in `projection_genes` does not match any rows of expression data")
        } else {
            message(
                sprintf("    Using supplied list of genes: Found %i/%i matches",
                    length(projection_genes), length(object@params$latentSpace$projectionGenes)
                    )
                )
        }
    }

    message()

    object@params$latentSpace$projectionGenes <- projection_genes


    return(object)
}


#' calculate signature scores
#'
#' For each signature-cell pair, compute a score that captures the level of
#' correspondence between the cell and the signature.
#' To estimate significance of these scores, a set of random gene signatures is
#' generated to create a null distribution
#'
#' @param object the VISION object
#' @param sigData a list of Signature objects for which to compute the scores
#' @param metaData a list of existing cell-signatures
#' @param sig_norm_method (optional) Method to apply to normalize the expression
#' matrix before calculating signature scores
#' @return the VISION object, with signature score slots populated
calcSignatureScores <- function(
    object, sigData = NULL,
    metaData = NULL, sig_norm_method = NULL) {

    message("Evaluating signature scores on cells...\n")

    ## override object parameters
    if (!is.null(sigData)) object@sigData <- sigData
    if (!is.null(metaData)) object@metaData <- metaData
    if (!is.null(sig_norm_method)) object@params$signatures$sigNormMethod <- sig_norm_method

    if (length(object@sigData) == 0) {
        sigScores <- matrix(nrow = ncol(object@exprData), ncol = 0,
                            dimnames = list(colnames(object@exprData), NULL)
                            )
        object@SigScores <- sigScores
        return(object)
    }

    normExpr <- getNormalizedCopySparse(
        object@exprData,
        object@params$signatures$sigNormMethod
    )

    sigScores <- batchSigEvalNorm(object@sigData, normExpr)

    object@SigScores <- sigScores

    return(object)
}


#' calculate gene-signature importance
#'
#' For each signature, the contribution of each gene to the signature score
#' is evaluated by calculating the covariance between signature scores and expression
#' The correlation of genes with a negative sign in the signature are inverted.
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colSds
#' @importFrom matrixStats rowSds
#'
#' @param object the VISION object
#' @return the VISION object, with SigGeneImportance slot populated
evalSigGeneImportance <- function(object){

    message("Evaluating signature-gene importance...\n")

    sigScores <- object@SigScores

    if (length(object@sigData) == 0) {
        object@SigGeneImportance <- list()
        return(object)
    }

    if (length(sigScores) <= 1){
        stop(
            sprintf("Signature scores have not yet been computed.  `calcSignatureScores` must be run before running `evalSigGeneImportance`")
            )
    }

    normExpr <- getNormalizedCopy(
        object@exprData,
        object@params$signatures$sigNormMethod)

    # Center each column of sigScores first

    mu <- colMeans(sigScores)

    sigScores <- t(sigScores)
    sigScores <- (sigScores - mu)
    sigScores <- t(sigScores)

    # Center each row of normExpr
    mu <- rowMeans(normExpr)

    normExpr <- (normExpr - mu)

    # Compute Covariances
    sigData <- object@sigData

    sigGene <- function(signame) {
        sigdata <- sigData[[signame]]

        genes <- sigdata@sigDict

        sigvals <- sigScores[, signame]

        geneIndices <- match(names(genes), rownames(normExpr))

        corr <- sigGeneInner(sigvals, normExpr, geneIndices)

        names(corr) <- names(genes)

        corr <- corr * genes

        return(corr)
    }

    sigs <- colnames(sigScores)
    res <- pbmclapply(setNames(sigs, sigs), sigGene)

    object@SigGeneImportance <- res


    return(object)
}


#' calculate gene-signature importance
#'
#' For each signature, the contribution of each gene to the signature score
#' is evaluated by calculating the covariance between signature scores and expression
#' The correlation of genes with a negative sign in the signature are inverted.
#'
#' This version is made to avoid inflating sparse matrices
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colSds
#' @importFrom matrixStats rowSds
#' @importFrom Matrix Matrix
#' @importFrom Matrix Diagonal
#'
#' @param object the VISION object
#' @return the VISION object, with SigGeneImportance slot populated
evalSigGeneImportanceSparse <- function(object){

    message("Evaluating signature-gene importance...\n")

    sigScores <- object@SigScores

    if (length(object@sigData) == 0) {
        object@SigGeneImportance <- list()
        return(object)
    }

    if (length(sigScores) <= 1){
        stop(
            sprintf("Signature scores have not yet been computed.  `calcSignatureScores` must be run before running `evalSigGeneImportance`")
            )
    }

    normExpr <- getNormalizedCopySparse(
        object@exprData,
        object@params$signatures$sigNormMethod)

    # Center each column of sigScores first

    mu <- colMeans(sigScores)

    sigScores <- t(sigScores)
    sigScores <- (sigScores - mu)
    sigScores <- t(sigScores)

    # Precompute some matrices we'll need later
    NGenes <- nrow(normExpr@data)
    NCells <- ncol(normExpr@data)
    Cog <- Matrix(1, ncol = 1, nrow = NGenes)
    Coc <- Matrix(normExpr@colOffsets, nrow = 1)
    Cs <- Diagonal(x = normExpr@colScaleFactors)
    Roc <- Matrix(1, nrow = 1, ncol = NCells)
    C1 <- t(Roc)

    RM <- normExpr@data %*% (Cs %*% C1) + Cog %*% (Coc %*% (Cs %*% C1))
    RM <- RM / NCells
    RM <- RM[, 1]

    # Compute Covariances
    sigData <- object@sigData

    sigGene <- function(signame) {
        sigdata <- sigData[[signame]]

        genes <- sigdata@sigDict

        S <- sigScores[, signame, drop = F]
        geneIndices <- rownames(normExpr@data) %in% names(genes)

        E <- normExpr@data[geneIndices, , drop = FALSE]

        Rog <- Matrix(RM[geneIndices], ncol = 1)
        Cog <- Matrix(1, ncol = 1, nrow = length(genes))

        geneCov <- E %*% Cs %*% S - Rog %*% (Roc %*% S) + Cog %*% (Coc %*% (Cs %*% S))
        geneCov <- geneCov / (NCells - 1)
        geneCov <- geneCov[, 1]
        geneCov <- geneCov[names(genes)]

        geneCov <- geneCov * genes # invert sign for negative genes

        return(geneCov)
    }

    sigs <- colnames(sigScores)
    res <- pbmclapply(setNames(sigs, sigs), sigGene)

    object@SigGeneImportance <- res

    return(object)
}


#' Computes the latent space of the expression matrix using PCA
#'
#' @param object the VISION object for which compute the latent space
#' @param projection_genes character vector of gene names to use for projections
#' @param perm_wPCA If TRUE, apply permutation wPCA to determine significant
#' number of components. Default is FALSE.
#' @return the VISION with latentSpace populated
computeLatentSpace <- function(object, projection_genes = NULL,
                               perm_wPCA = NULL) {

    message("Computing a latent space for expression data...\n")

    if (!is.null(projection_genes)) object@params$latentSpace$projectionGenes <- projection_genes
    if (!is.null(perm_wPCA)) object@params$latentSpace$permPCA <- perm_wPCA

    expr <- object@exprData
    projection_genes <- object@params$latentSpace$projectionGenes
    perm_wPCA <- object@params$latentSpace$permPCA

    if (!is.null(projection_genes)) {
        exprData <- expr[projection_genes, , drop = FALSE]
    } else {
        exprData <- expr
    }

    exprData <- matLog2(exprData)

    if (perm_wPCA) {
        res <- applyPermutationWPCA(exprData, components = 30)
        pca_res <- res[[1]]
    } else {
        res <- applyPCA(exprData, maxComponents = 30)
        pca_res <- res[[1]]
    }

    object@LatentSpace <- t(pca_res)
    colnames(object@LatentSpace) <- paste0("PC ", seq_len(ncol(object@LatentSpace)))
    object@params$latentSpace$name <- "PCA"
    return(object)
}

#' generate projections
#'
#' Generates 2-dimensional representations of the expression matrix
#' Populates the 'Projections' slot on the VISION object
#'
#'
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
generateProjections <- function(object) {
  message("Projecting data into 2 dimensions...")

  # Some projection methods operate on the full expression matrix
  # If using one of these, we need to compute 'projection_genes'
  projection_methods <- object@params$projectionMethods
  if ("ICA" %in% projection_methods || "RBFPCA" %in% projection_methods) {
      object <- filterData(object)
  }

  projections <- generateProjectionsInner(object@exprData,
                                     object@LatentSpace,
                                     projection_genes = object@params$latentSpace$projectionGenes,
                                     projection_methods = object@params$projectionMethods,
                                     K = object@params$numNeighbors)

  # Add already-input projections
  for (proj in names(object@Projections)){
      projections[[proj]] <- object@Projections[[proj]]
  }

  # Make sure all projections have column names
  n <- names(projections)
  projections <- lapply(setNames(n, n), function(pname){
      proj <- projections[[pname]]
      if (is.null(colnames(proj))){
          colnames(proj) <- paste0(pname, "-", seq_len(ncol(proj)))
      }
      return(proj)
  })

  object@Projections <- projections

  message("")

  return(object)
}

#' Compute local correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
analyzeLocalCorrelations <- function(object) {

  signatureBackground <- generatePermutationNull(
      object@exprData, object@sigData, num = 3000
  )

  normExpr <- getNormalizedCopySparse(
      object@exprData,
      object@params$signatures$sigNormMethod)

  message("Evaluating local consistency of signatures in latent space...\n")

  weights <- computeKNNWeights(object@LatentSpace, object@params$numNeighbors)

  sigConsistencyScores <- sigConsistencyScores(
                                weights,
                                object@SigScores,
                                object@metaData,
                                signatureBackground,
                                normExpr)

  if (hasFeatureBarcodeData(object)) {
      fbcs <- fbConsistencyScores(weights, object@featureBarcodeData)
  } else {
      fbcs <- NULL
  }

  message("Clustering signatures...\n")
  sigClusters <- clusterSignatures(object@SigScores,
                                   object@metaData,
                                   sigConsistencyScores,
                                   clusterMeta = object@params$micropooling$pool)

  metaConsistencyScores <- sigConsistencyScores[
      colnames(object@metaData), , drop = FALSE
  ]

  sigConsistencyScores <- sigConsistencyScores[
      colnames(object@SigScores), , drop = FALSE
  ]


  LocalAutocorrelation <- list(
      "Signatures" = sigConsistencyScores,
      "Meta" = metaConsistencyScores,
      "Features" = fbcs,
      "Clusters" = sigClusters
  )

  object@LocalAutocorrelation <- LocalAutocorrelation

  return(object)
}


#' Compute trajectory correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
analyzeTrajectoryCorrelations <- function(object) {

  signatureBackground <- generatePermutationNull(
      object@exprData, object@sigData, num = 3000
  )

  normExpr <- getNormalizedCopySparse(
      object@exprData,
      object@params$signatures$sigNormMethod)

  message("Evaluating local consistency of signatures within trajectory model...\n")

  weights <- computeKNNWeights(object@LatentTrajectory, object@params$numNeighbors)

  sigVTreeProj <- sigConsistencyScores(weights,
                                       object@SigScores,
                                       object@metaData,
                                       signatureBackground,
                                       normExpr)

  if (hasFeatureBarcodeData(object)) {
      fbcs <- fbConsistencyScores(weights, object@featureBarcodeData)
  } else {
      fbcs <- NULL
  }

  message("Clustering signatures...\n")
  sigTreeClusters <- clusterSignatures(object@SigScores,
                                       object@metaData,
                                       sigVTreeProj,
                                       clusterMeta = object@params$micropooling$pool)

  metaConsistencyScores <- sigVTreeProj[
      colnames(object@metaData), , drop = FALSE
  ]

  sigConsistencyScores <- sigVTreeProj[
      colnames(object@SigScores), , drop = FALSE
  ]

  TrajectoryAutocorrelation <- list(
      "Signatures" = sigConsistencyScores,
      "Meta" = metaConsistencyScores,
      "Features" = fbcs,
      "Clusters" = sigClusters
  )

  object@TrajectoryAutocorrelation <- TrajectoryAutocorrelation

  return(object)
}


#' Compute Ranksums Test, for all factor meta data.  One level vs all others
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colRanks
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
clusterSigScores <- function(object) {

    message("Computing differential signature tests...\n")

    sigScores <- object@SigScores
    metaData <- object@metaData

    metaData <- metaData[rownames(sigScores), , drop = FALSE]

    # Determine which metaData we can run on
    # Must be a factor with at least 20 levels
    clusterMeta <- vapply(colnames(metaData), function(x) {
            scores <- metaData[[x]]
            if (!is.factor(scores)){
                return("")
            }
            if (length(levels(scores)) > 50){
                return("")
            }
            if (length(unique(scores)) == 1){
                return("")
            }
            return(x)
        }, FUN.VALUE = "")
    clusterMeta <- clusterMeta[clusterMeta != ""]

    ClusterComparisons <- list()

    # Comparisons for Signatures
    if (ncol(sigScores) > 0){
        sigScoreRanks <- colRanks(sigScores,
                                  preserveShape = TRUE,
                                  ties.method = "average")
        dimnames(sigScoreRanks) <- dimnames(sigScores)
    } else {
        sigScoreRanks <- sigScores
    }

    out <- pbmclapply(clusterMeta, function(variable){
        values <- metaData[[variable]]
        var_levels <- levels(values)

        result <- lapply(var_levels, function(var_level){
            cluster_ii <- which(values == var_level)

            r1 <- matrix_wilcox(sigScoreRanks, cluster_ii,
                                check_na = FALSE, check_ties = FALSE)

            pval <- r1$pval
            stat <- r1$stat
            fdr <- p.adjust(pval, method = "BH")
            out <- data.frame(
                stat = stat, pValue = pval, FDR = fdr
            )
            return(out)
        })

        names(result) <- var_levels
        result <- result[order(var_levels)]

        return(result)
    }, mc.cores = 1)

    ClusterComparisons[["Signatures"]] <- out

    # Comparisons for Meta data
    # Split meta into numeric and factor
    numericMeta <- vapply(seq_len(ncol(metaData)),
                          function(i) is.numeric(metaData[[i]]),
                          FUN.VALUE = TRUE)
    numericMeta <- metaData[, numericMeta, drop = F]
    numericMeta <- as.matrix(numericMeta)

    factorMeta <- vapply(seq_len(ncol(metaData)),
                          function(i) is.factor(metaData[[i]]),
                          FUN.VALUE = TRUE)
    factorMeta <- metaData[, factorMeta, drop = F]

    if (ncol(numericMeta) > 0){
        numericMetaRanks <- colRanks(numericMeta,
                                     preserveShape = TRUE,
                                     ties.method = "average")
        dimnames(numericMetaRanks) <- dimnames(numericMeta)
    } else {
        numericMetaRanks <- numericMeta
    }

    out <- pbmclapply(clusterMeta, function(variable){
        values <- metaData[[variable]]
        var_levels <- levels(values)

        result <- lapply(var_levels, function(var_level){
            cluster_ii <- which(values == var_level)

            r2 <- matrix_wilcox(numericMetaRanks, cluster_ii,
                                check_na = TRUE, check_ties = TRUE)

            r3 <- matrix_chisq(factorMeta, cluster_ii)

            pval <- c(r2$pval, r3$pval)
            stat <- c(r2$stat, r3$stat)
            fdr <- p.adjust(pval, method = "BH")
            out <- data.frame(
                stat = stat, pValue = pval, FDR = fdr
            )
            return(out)
        })

        names(result) <- var_levels
        result <- result[order(var_levels)]

        return(result)
    }, mc.cores = 1)

    ClusterComparisons[["Meta"]] <- out

    # Comparisons for Features
    if (hasFeatureBarcodeData(object)){
        fbRanks <- colRanks(object@featureBarcodeData,
                            preserveShape = TRUE,
                            ties.method = "average")
        dimnames(fbRanks) <- dimnames(object@featureBarcodeData)
        out_fb <- pbmclapply(clusterMeta, function(variable){
            values <- metaData[[variable]]
            var_levels <- levels(values)

            result <- lapply(var_levels, function(var_level){
                cluster_ii <- which(values == var_level)

                rr <- matrix_wilcox(fbRanks, cluster_ii,
                                    check_na = TRUE, check_ties = TRUE)

                pval <- rr$pval
                stat <- rr$stat
                fdr <- p.adjust(pval, method = "BH")
                out <- data.frame(
                    stat = stat, pValue = pval, FDR = fdr
                )
                return(out)
            })

            names(result) <- var_levels
            result <- result[order(var_levels)]

            return(result)
        }, mc.cores = 1)

        ClusterComparisons[["Features"]] <- out_fb
    } else {
        ClusterComparisons[["Features"]] <- NULL
    }

    object@ClusterComparisons <- ClusterComparisons

    return(object)

}

#' Compute pearson correlation between signature scores and principle components
#'
#' Populations the PCAnnotatorData slot of the VISION object
#'
#' @importFrom pbmcapply pbmclapply
#' @param object the VISION object
#' @return pearsonCorr numeric matrix N_Signatures x N_PCs
calculatePearsonCorr <- function(object){

  message("Computing correlations between signatures and latent space components...\n")
  sigMatrix <- object@SigScores
  metaData <- object@metaData
  latentSpace <- object@LatentSpace

  ## combined gene signs and numeric meta variables

  numericMetaVars <- vapply(colnames(metaData),
                            function(x) is.numeric(metaData[[x]]),
                            FUN.VALUE = TRUE)

  numericMeta <- metaData[, numericMetaVars, drop = FALSE]
  numericMeta <- numericMeta[rownames(sigMatrix), , drop = FALSE]


  computedSigMatrix <- cbind(sigMatrix, numericMeta)

  pearsonCorr <- pbmclapply(seq_len(ncol(computedSigMatrix)), function(i) {
      ss <- computedSigMatrix[, i];

      ls_col_cor <- apply(latentSpace, 2, function(pc){
           suppressWarnings({
               pc_result <- cor.test(ss, pc)
           })
           if (is.na(pc_result$estimate)) {  # happens if std dev is 0 for a sig
               return(0)
           } else {
               return(pc_result$estimate)
           }
      })

      return(ls_col_cor)
  })

  pearsonCorr <- do.call(rbind, pearsonCorr)
  rownames(pearsonCorr) <- colnames(computedSigMatrix)
  colnames(pearsonCorr) <- colnames(latentSpace)

  if (hasFeatureBarcodeData(object)) {
      featureBarcodeData <- object@featureBarcodeData
      pearsonCorrFeatures <- pbmclapply(
          seq_len(ncol(featureBarcodeData)), function(i) {
          ss <- featureBarcodeData[, i];

          ls_col_cor <- apply(latentSpace, 2, function(pc){
               suppressWarnings({
                   pc_result <- cor.test(ss, pc)
               })
               if (is.na(pc_result$estimate)) {  # happens if std dev is 0 for a feature
                   return(0)
               } else {
                   return(pc_result$estimate)
               }
          })

          return(ls_col_cor)
      })

      pearsonCorrFeatures <- do.call(rbind, pearsonCorrFeatures)
      rownames(pearsonCorrFeatures) <- colnames(featureBarcodeData)
      colnames(pearsonCorrFeatures) <- colnames(latentSpace)
  } else {
      pearsonCorrFeatures <- NULL
  }


  pcaAnnotData <- PCAnnotatorData(
      pearsonCorr = pearsonCorr, pearsonCorrFeatures = pearsonCorrFeatures
  )
  object@PCAnnotatorData <- pcaAnnotData

  return(object)
}
