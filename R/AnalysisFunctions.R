#' Creates an initial clustering of the cells
#'
#' Results of this are stored as a new variabe in the object's metaData
#' and 'cluster_variable' is populated with its name
#'
#' @param object the VISION object for which to cluster the cells
#' @return the VISION object modifed as described above
clusterCells <- function(object) {

    message("Clustering cells...")

    if (sum(dim(object@latentSpace)) == 2) { # No latent Space

        exprData <- matLog2(object@exprData)

        filterInput <- object@projection_genes
        filterThreshold <- object@threshold

        gene_passes <- applyFilters(exprData, filterThreshold, filterInput)
        fexpr <- exprData[gene_passes, ]

        # Compute wcov using matrix operations to avoid
        # creating a large dense matrix

        N <- ncol(fexpr)
        wcov <- tcrossprod(fexpr) / N

        mu <- as.matrix(rowMeans(fexpr), ncol = 1)
        mumu <- tcrossprod(mu)
        wcov <- as.matrix(wcov - mumu)

        # SVD of wieghted correlation matrix
        ncomp <- min(ncol(fexpr), nrow(fexpr), 10)
        decomp <- rsvd::rsvd(wcov, k = ncomp)
        evec <- t(decomp$u)

        # Project down using computed eigenvectors
        res <- (evec %*% fexpr) - as.vector(evec %*% mu)
        res <- as.matrix(res)
        res <- t(res) # avoid transposing many times below
    } else {
        res <- object@latentSpace
    }

    n_workers <- getWorkerCount()

    kn <- ball_tree_knn(res,
                        min(round(sqrt(nrow(res))), 30),
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
    object@cluster_variable <- cluster_variable

    return(object)

}


#' create micro-clusters that reduce noise and complexity while maintaining
#' the overall signal in the data
#' @param object the VISION object for which to cluster the cells
#' @param cellsPerPartition the minimum number of cells to put into a cluster
#' @return the VISION with pooled cells
poolCells <- function(object,
                      cellsPerPartition=object@cellsPerPartition) {
    object@cellsPerPartition <- cellsPerPartition

    message(paste(
      "Performing micro-pooling on",
      ncol(object@exprData),
      "cells with a target pool size of",
      object@cellsPerPartition
    ))

    if (object@cluster_variable != "") {
        preserve_clusters <- object@metaData[[object@cluster_variable]]
        names(preserve_clusters) <- rownames(object@metaData)
    } else {
        preserve_clusters <- NULL
    }

    pools <- applyMicroClustering(object@exprData,
                                          cellsPerPartition = object@cellsPerPartition,
                                          filterInput = object@projection_genes,
                                          filterThreshold = object@threshold,
                                          preserve_clusters = preserve_clusters,
                                          latentSpace = object@latentSpace)

    object@pools <- pools

    pooled_cells <- createPoolsBatch(object@pools, object@exprData)
    object@exprData <- pooled_cells

    pooled_unnorm <- createPoolsBatch(object@pools, object@unnormalizedData)
    object@unnormalizedData <- pooled_unnorm

    if (!all(dim(object@latentSpace) == c(1, 1))) {
        pooled_latent <- t(createPoolsBatch(object@pools, t(object@latentSpace)))
        object@latentSpace <- pooled_latent
    }

    poolMeta <- createPooledMetaData(object@metaData, object@pools)
    object@metaData <- poolMeta

    if (length(object@inputProjections) > 0){

        newInputProjections <- lapply(
            object@inputProjections,
            function(proj) {
                new_coords <- t(createPoolsBatch(object@pools, t(proj)))
                return(new_coords)
            })

        names(newInputProjections) <- names(object@inputProjections)

        object@inputProjections <- newInputProjections

    }

    return(object)
}

#' filter data accourding to the provided filters
#' @param object the VISION object
#' @param threshold threshold to apply for the threshold filter
#' @param projection_genes either a list of genes or a method to select genes
#' @return the VISION object, populated with filtered data
filterData <- function(object,
                       threshold=object@threshold,
                       projection_genes=object@projection_genes) {

    object@projection_genes <- projection_genes
    object@threshold <- threshold

    message("Determining Projection Genes...")

    exprData <- matLog2(object@exprData)

    object@projection_genes <- applyFilters(
                exprData,
                object@threshold,
                object@projection_genes)

    return(object)
}

#' Calculate weights based on the false negative rates, computed using the
#' provided housekeeping genes
#' @param object the VISION object
#' @param nomodel (optional) if TRUE, no fnr curve calculated and all weights
#' equal to 1. Else FNR and weights calculated.
#' @return the VISION object with populated weights slot
calcWeights <- function(object,
                        nomodel=object@nomodel) {
    object@nomodel <- nomodel
    clustered <- (ncol(object@exprData) > 15000 || object@pool)
    if (!clustered && !object@nomodel) {
        message("Computing weights from False Negative Rate curves...")
        falseneg_out <- createFalseNegativeMap(object@exprData,
                                               object@housekeepingData)

        object@weights <- computeWeights(falseneg_out[[1]], falseneg_out[[2]],
                                         object@exprData)

        rownames(object@weights) <- rownames(object@exprData)
        colnames(object@weights) <- colnames(object@exprData)

    }

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
#' @param sig_score_method the scoring method to use
#' @return the VISION object, with signature score slots populated
calcSignatureScores <- function(object,
                                sigData=object@sigData,
                                metaData=object@metaData,
                                sig_norm_method=object@sig_norm_method,
                                sig_score_method=object@sig_score_method) {

    message("Evaluating Signature Scores on Cells...")

    ## override object parameters
    if(!is.null(sigData)) object@sigData <- sigData
    if(!is.null(metaData)) object@metaData <- metaData
    object@sig_norm_method <- sig_norm_method
    object@sig_score_method <- sig_score_method

    normExpr <- getNormalizedCopy(object@exprData, object@sig_norm_method)

    sigScores <- batchSigEval(object@sigData, object@sig_score_method,
                              normExpr, object@weights)

    object@sigScores <- sigScores

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

    message("Computing a latent space for expression data...")

    if (!is.null(projection_genes)) object@projection_genes <- projection_genes
    if (!is.null(perm_wPCA)) object@perm_wPCA <- perm_wPCA

    expr <- object@exprData
    weights <- object@weights
    projection_genes <- object@projection_genes
    perm_wPCA <- object@perm_wPCA

    if (!is.null(projection_genes)) {
        exprData <- expr[projection_genes, ]
    } else {
        exprData <- expr
    }

    if (all(dim(weights) == c(1, 1))){
        weights <- matrix(1, nrow = nrow(exprData), ncol = ncol(exprData),
                          dimnames = list(
                                          rownames(exprData),
                                          colnames(exprData)
                                         )
                         )
    }

    exprData <- matLog2(exprData)

    if (perm_wPCA) {
        res <- applyPermutationWPCA(exprData, weights, components = 30)
        pca_res <- res[[1]]
    } else {
        res <- applyWeightedPCA(exprData, weights, maxComponents = 30)
        pca_res <- res[[1]]
    }

    object@latentSpace <- t(pca_res)
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

  projections <- generateProjectionsInner(object@exprData,
                                     object@latentSpace,
                                     projection_genes = object@projection_genes,
                                     projection_methods = object@projection_methods)

  # Add inputProjections
  for (proj in names(object@inputProjections)){
      projections[[proj]] <- object@inputProjections[[proj]]
  }

  object@Projections <- projections

  return(object)
}

#' Compute spatial correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @param signatureBackground as returned by `calculateSignatureBackground`
#' @return the VISION object with values set for the analysis results
analyzeSpatialCorrelations <- function(object, signatureBackground = NULL) {

  if (is.null(signatureBackground)){
      message("Computing background distribution for signature scores...")
      signatureBackground <- calculateSignatureBackground(object, num = 3000)
  }

  message("Evaluating spatial consistency of signatures...")

  sigConsistencyScores <- sigConsistencyScores(
                                object@latentSpace,
                                object@sigScores,
                                object@metaData,
                                signatureBackground)

  message("Clustering Signatures...")
  sigClusters <- clusterSignatures(object@sigScores,
                                   object@metaData,
                                   sigConsistencyScores$emp_pVals,
                                   sigConsistencyScores$sigProjMatrix,
                                   clusterMeta = object@pool)

  sigConsistencyScoresData <- ProjectionData(
                             sigProjMatrix = sigConsistencyScores$sigProjMatrix,
                             pMatrix = sigConsistencyScores$pVals,
                             sigClusters = sigClusters,
                             emp_pMatrix = sigConsistencyScores$emp_pVals)

  object@SigConsistencyScores <- sigConsistencyScoresData

  return(object)
}


#' Compute trajectory correlations for all signatures
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured succesfully by the projections.
#' @param object the VISION object
#' @param signatureBackground as returned by `calculateSignatureBackground`
#' @return the VISION object with values set for the analysis results
analyzeTrajectoryCorrelations <- function(object, signatureBackground = NULL) {

  if (is.null(signatureBackground)){
      message("Computing background distribution for signature scores...")
      signatureBackground <- calculateSignatureBackground(object, num = 3000)
  }

  message("Computing significance of signatures...")
  sigVTreeProj <- sigConsistencyScores(object@latentTrajectory,
                                       object@sigScores,
                                       object@metaData,
                                       signatureBackground)

  message("Clustering Signatures...")
  sigTreeClusters <- clusterSignatures(object@sigScores,
                                       object@metaData,
                                       sigVTreeProj$pVals,
                                       sigVTreeProj$sigProjMatrix,
                                       clusterMeta = object@pool)

  TrajectoryConsistencyScores <- ProjectionData(
      sigProjMatrix = sigVTreeProj$sigProjMatrix,
      pMatrix = sigVTreeProj$pVals,
      sigClusters = sigTreeClusters,
      emp_pMatrix = sigVTreeProj$emp_pVals)

  object@TrajectoryConsistencyScores <- TrajectoryConsistencyScores

  return(object)
}


#' Compute Ranksums Test, for all factor meta data.  One level vs all others
#'
#' @importFrom parallel mclapply
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
clusterSigScores <- function(object) {

    message("Computing differential signature tests...")

    sigScores <- object@sigScores
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

    clusterPVals <- function(sigScores, metaData, variable) {

        values <- metaData[[variable]]
        var_levels <- levels(values)

        result <- lapply(var_levels, function(var_level){

            cluster_ii <- which(values == var_level)
            not_cluster_ii <- which(values != var_level)

            # Process the gene signatures
            pvals <- mclapply(colnames(sigScores), function(sig){

                if(length(cluster_ii) == 0 || length(not_cluster_ii) == 0){
                    return(list(pval = 1.0, stat = 0))
                }

                suppressWarnings({
                    out <- wilcox.test(sigScores[cluster_ii, sig],
                                   sigScores[not_cluster_ii, sig],
                                   alternative = "two.sided", exact = FALSE)
                })
                pval <- out$p.value
                AUC <- out$statistic / length(cluster_ii) / length(not_cluster_ii)
                stat <- (AUC - .5) * 2  # translates (0.5, 1) to (0, 1)
                return(list(pval = pval, stat = stat))
            }, mc.cores = 10)
            names(pvals) <- colnames(sigScores)

            # Process the metaData variables
            meta_pvals <- lapply(colnames(metaData), function(sig){

                if (length(cluster_ii) == 0 || length(not_cluster_ii) == 0){
                    return(list(pval = 1.0, stat = 0))
                }

                if (is.numeric(metaData[, sig])) {

                    x <- metaData[cluster_ii, sig]
                    y <- metaData[not_cluster_ii, sig]

                    if (all(is.na(x)) || all(is.na(y))){
                        pval <- 1.0
                        stat <- 0
                    } else {

                        suppressWarnings({
                            out <- wilcox.test(x, y,
                                   alternative = "two.sided", exact = FALSE)
                        })

                        pval <- out$p.value
                        AUC <- out$statistic / length(cluster_ii) / length(not_cluster_ii)
                        stat <- (AUC - .5) * 2  # translates (0, 1) to (-1, 1)

                    }
                } else if (is.factor(metaData[, sig])) {
                    sigvals <- metaData[, sig, drop = F]
                    sigvals[, 2] <- 0
                    sigvals[cluster_ii, 2] <- 1
                    sigvals[, 2] <- as.factor(sigvals[, 2])
                    M <- table(sigvals)
                    if (nrow(M) > 1){
                        suppressWarnings(
                            out <- chisq.test(M)
                        )
                        pval <- out$p.value
                        n <- sum(M)
                        V <- sqrt(out$statistic / n /
                            min(nrow(M) - 1, ncol(M) - 1)
                        )
                        stat <- V # Cramer's V
                    } else {
                        pval <- 1.0
                        stat <- 0
                    }
                }
                return(list(pval = pval, stat = stat))
            })
            names(meta_pvals) <- colnames(metaData)

            all <- c(pvals, meta_pvals)
            all_pvals <- vapply(all, function(x) x$pval, FUN.VALUE = 0.0)
            all_zscores <- vapply(all, function(x) x$stat, FUN.VALUE = 0.0)

            return(list(pvals = all_pvals, zscores = all_zscores))
        })
        names(result) <- var_levels
        zscores <- do.call(cbind,
                         lapply(result, function(x) x$zscores))
        pvals <- do.call(cbind,
                         lapply(result, function(x) x$pvals))

        pvals_adj <- apply(pvals, MARGIN = 2, FUN = p.adjust, method = "BH")
        pvals_adj[pvals_adj < 1e-300] <- 1e-300
        pvals_adj <- log10(pvals_adj)
        return(list(pvals = pvals_adj, zscores = zscores))
    }

    pvalsAndZscores <- lapply(clusterMeta, function(variable){
        result <- clusterPVals(sigScores, metaData, variable)
        return(result)
    })

    object@ClusterSigScores <- pvalsAndZscores
    return(object)

}

#' Compute pearson correlation between signature scores and principle components
#'
#' Populations the PCAnnotatorData slot of the VISION object
#'
#' @param object the VISION object
#' @return pearsonCorr numeric matrix N_Signatures x N_PCs
calculatePearsonCorr <- function(object){

  message("Computing Correlations between Signatures and Expression PCs...")
  sigMatrix <- object@sigScores
  metaData <- object@metaData
  latentSpace <- object@latentSpace

  ## combined gene signs and numeric meta variables

  numericMetaVars <- vapply(colnames(metaData),
                            function(x) is.numeric(metaData[[x]]),
                            FUN.VALUE = TRUE)

  numericMeta <- metaData[, numericMetaVars, drop = FALSE]
  numericMeta <- numericMeta[rownames(sigMatrix), , drop = FALSE]


  computedSigMatrix <- cbind(sigMatrix, numericMeta)

  pearsonCorr <- lapply(1:ncol(computedSigMatrix), function(i) {
    lapply(1:ncol(latentSpace), function(j) {
               ss <- computedSigMatrix[, i];
               pc <- latentSpace[, j];
               suppressWarnings({
                   pc_result <- cor.test(ss, pc)
               })
               if (is.na(pc_result$estimate)) {  # happens is std dev is 0 for a sig
                   return(0)
               } else {
                   return(pc_result$estimate)
               }
    })
  })

  pearsonCorr <- matrix(unlist(lapply(pearsonCorr, unlist)),
                        nrow=ncol(computedSigMatrix),
                        ncol=ncol(latentSpace), byrow=TRUE)

  rownames(pearsonCorr) <- colnames(computedSigMatrix)
  colnames(pearsonCorr) <- colnames(latentSpace)

  pcaAnnotData <- PCAnnotatorData(pearsonCorr = pearsonCorr)
  object@PCAnnotatorData <- pcaAnnotData

  return(object)
}

convertToDense <- function(object) {

    object@exprData <- as.matrix(object@exprData)

    return(object)
}
