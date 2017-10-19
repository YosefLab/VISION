# require(geometry)
# require(pROC)

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
            .Object <- new("Signature")
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
#'
#' @param N_SAMPLES Number of samples to generate the distribution for
#' @param NUM_REPLICATES Number of replicates to generate for the background distribution
#' @return Random matrix with dimensions N_SAMPLES x NUM_REPLICATES with row ordered in ascending order
getBGDist <- function(N_SAMPLES, NUM_REPLICATES) {

    if (nrow(BG_DIST) != N_SAMPLES || ncol(BG_DIST) != NUM_REPLICATES) {
    set.seed(RANDOM_SEED)
    BG_DIST <- matrix(stats::rnorm(N_SAMPLES*NUM_REPLICATES),
                        nrow=N_SAMPLES, ncol=NUM_REPLICATES)
    BG_DIST <- apply(BG_DIST, 2, order)
    }

    return(BG_DIST)

}

#' Evaluates the significance of each signature vs. each projection.
#' @importFrom Matrix crossprod
#' @importFrom mclust densityMclust
#' @importFrom entropy entropy.plugin
#' @importFrom pROC multiclass.roc
#' @param projections Maps projections to their spatial coordinates for each
#' sample
#' @param sigScoresData List of SignatureScores Object, mapping signature
#' names to their value at each coordinate
#' @param randomSigData List of SignatureScores Object, mapping randomly
#' generated signatures to scores to be compared with the real signature scores.
#' @param BPPARAM the parallelization backend to use
#' @return list:
#' \itemize{
#'     \item sigNames: labels for the signatures (rows in output matrices)
#'     \item projNames: labels for the projections (columns in output matrices)
#'     \item sigProbMatrix: the matrix of signature-projection consistency scores
#'     \item pVals: pvalues for the scores
#' }
sigsVsProjections <- function(projections, sigScoresData,
                                randomSigData, BPPARAM=bpparam()) {

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

    factorSigProjMatrix <- matrix(0L, nrow=N_SIGNATURE_FACTORS,
                                ncol=N_PROJECTIONS)
    factorSigProjMatrix_P <- matrix(0L, nrow=N_SIGNATURE_FACTORS,
                                    ncol=N_PROJECTIONS)

    pnumSigProjMatrix <- matrix(0L, nrow=N_SIGNATURE_PNUM, ncol=N_PROJECTIONS)
    pnumSigProjMatrix_P <- matrix(0L, nrow=N_SIGNATURE_PNUM, ncol=N_PROJECTIONS)

    # Build a matrix of all signatures
    sigScoreMatrix <- matrix(unlist(bplapply(spRows, function(sig) {
    rank(sig@scores, ties.method="average")
    },
    BPPARAM=BPPARAM)), nrow=N_SAMPLES, ncol=length(spRows))


    randomSigScoreMatrix <- matrix(unlist(bplapply(randomSigData, function(rsig) {
    rank(rsig@scores, ties.method="average")
    },
    BPPARAM=BPPARAM)), nrow=N_SAMPLES, ncol=length(randomSigData))


    ### build one hot matrix for factors
    factorSigs <- list()
    for (s in precomputedFactor) {
    fValues <- s@scores
    fLevels <- levels(fValues)
    factorFreq <- matrix(0L, ncol=length(fLevels))
    factorMatrix <- matrix(0L, nrow=N_SAMPLES, ncol=length(fLevels))

    factList <- lapply(fLevels, function(fval) {
                factorMatrixRow <- matrix(0L, nrow=N_SAMPLES, ncol=1)
                equal_ii <- which(fValues == fval)
                factorMatrixRow[equal_ii] <- 1
                return(list(length(equal_ii) / length(fValues), factorMatrixRow))
            })
    factorFreq <- lapply(factList, function(x) return(x[[1]]))
    factorMatrix <- matrix(unlist(lapply(factList, function(x) return(x[[2]]))),
                            nrow=N_SAMPLES, ncol=length(fLevels))

    factorSigs[[s@name]] <- list(fLevels, factorFreq, factorMatrix)
    }

    message("Evaluating signatures against projections...")

    i <- 1
    projnames <- names(projections)
    for (proj in projections) {
        weights <- computeKNNWeights(proj, K=round(sqrt(NCOL(proj@pData))), BPPARAM)
        neighborhoodPrediction <- Matrix::crossprod(weights, sigScoreMatrix)

        ## Neighborhood dissimilatory score = |actual - predicted|
        dissimilarity <- abs(sigScoreMatrix - neighborhoodPrediction)
        medDissimilarity <- as.matrix(apply(dissimilarity, 2, stats::median))

    # Calculate scores for random signatures
    randomNeighborhoodPrediction <- Matrix::crossprod(weights,
                                                        randomSigScoreMatrix)
    randomDissimilarity <- abs(randomSigScoreMatrix -
                                    randomNeighborhoodPrediction)
    randomMedDissimilarity <- as.matrix(apply(randomDissimilarity, 2,
                                                stats::median))

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

        bgStat <- matrix(unlist(bplapply(names(backgrounds), function(x) {
                mu_x <- mean(backgrounds[[x]])
                std_x <- biasedVectorSD(as.matrix(backgrounds[[x]]))
                return(list(as.numeric(x), mu_x, std_x))
            }, BPPARAM=BPPARAM)), nrow=length(names(backgrounds)),
                ncol=3, byrow=TRUE)


        mu <- matrix(unlist(bplapply(spRows, function(x) {
                numG <- x@numGenes
                row_i <- which.min(abs(numG - bgStat[,1]))
                return(bgStat[row_i, 2])
            }, BPPARAM=BPPARAM)), nrow=nrow(medDissimilarity),
                ncol=ncol(medDissimilarity))

        sigma <- matrix(unlist(bplapply(spRows, function(x) {
                numG <- x@numGenes
                row_i <- which.min(abs(numG - bgStat[,1]))
                return(bgStat[row_i,3])
            }, BPPARAM=BPPARAM)), nrow=nrow(medDissimilarity),
                ncol=ncol(medDissimilarity))


        #Create CDF function for medDissmilarityPrime and apply CDF function to
        #medDissimilarityPrime pointwise
        pValues <- stats::pnorm( ((medDissimilarity - mu) / sigma))

        sigProjMatrix[,i] <- 1 - (medDissimilarity / N_SAMPLES)
        sigProjMatrix_P[,i] <- pValues

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

                krTest <- stats::kruskal.test(krList)

                return(list(a, krTest$p.value))
            })

        factorSigProjMatrix[,i] <- unlist(lapply(factorSigProjList,
                                                function(x) x[[1]]))
        factorSigProjMatrix_P[,i] <- unlist(lapply(factorSigProjList,
                                                    function(x) x[[2]]))

    i <- i+1
    }

    # Concatenate Factor Sig-proj entries back in
    sigProjMatrix <- rbind(sigProjMatrix, factorSigProjMatrix, pnumSigProjMatrix)
    sigProjMatrix_P <- rbind(sigProjMatrix_P,
                            factorSigProjMatrix_P,
                            pnumSigProjMatrix_P)
    spRowLabels <- c(spRowLabels, spRowLabelsFactors, spRowLabelsPNum)

    original_shape <- dim(sigProjMatrix_P)
    sigProjMatrix_P <- matrix(stats::p.adjust(sigProjMatrix_P, method="BH"),
                            nrow=nrow(sigProjMatrix_P),
                            ncol=ncol(sigProjMatrix_P))
    sigProjMatrix_P[sigProjMatrix_P == 0] <- 10^(-300)

    sigProjMatrix_P <- as.matrix(log10(sigProjMatrix_P))

    colnames(sigProjMatrix_P) <- projnames
    rownames(sigProjMatrix_P) <- spRowLabels

    colnames(sigProjMatrix) <- projnames
    rownames(sigProjMatrix) <- spRowLabels

    return(list(sigNames = spRowLabels, projNames = spColLabels,
                sigProjMatrix = sigProjMatrix, pVals = sigProjMatrix_P))
}

#' Clusters signatures according to the rank sum
#' @importFrom mclust densityMclust
#' @param sigList List of signatures
#' @param sigMatrix Matrix of signatures scores, NUM_SIGNATURES x NUM_SAMPLES
#' @param pvals the corresponding P-values for each score,
#' NUM_SIGNATURES x NUM_SAMPLES
#' @param k Number of clusters to generate
#' @return a list:
#' \itemize{
#'     \item Computed: a list of clusters of computed signatures
#'     \item Precomputed: a list of clusters of precomputed signatures
#' }
clusterSignatures <- function(sigList, sigMatrix, pvals, k=10) {

    precomputed <- lapply(sigList, function(x) x@isPrecomputed)

    significant <- apply(pvals, 1, function(x) min(x) < -1.3)

    signif_computed <- significant[names(which(precomputed == FALSE))]
    signif_precomp <- significant[names(which(precomputed == TRUE))]

    keep_computed = names(which(signif_computed == TRUE))

    # Cluster computed signatures and precomputed signatures separately
    computedSigsToCluster <- names(precomputed[which(precomputed==FALSE)])
    computedSigMatrix <- sigMatrix[computedSigsToCluster,,drop=FALSE]

    computedSigMatrix <- computedSigMatrix[keep_computed,,drop=FALSE]

    compcls <- list()
    maxcls <- 1
    if (nrow(computedSigMatrix) > 1) {
        r <- as.matrix(t(apply(computedSigMatrix, 1,
                                function(x) rank(x, ties.method="average"))))

    compkm <- densityMclust(r)
    compcls <- as.list(compkm$classification)
    compcls <- compcls[order(unlist(compcls), decreasing=FALSE)]

    maxcls <- max(unlist(compcls))
    }

    compcls[names(which(significant == FALSE))] <- maxcls + 1

    # Don't actually cluster Precomputed Signatures -- just return in a list.
    precompcls <- list()
    if (length(which(precomputed==TRUE)) > 0) {
        precomputedSigsToCluster <- names(precomputed[which(precomputed==TRUE)])
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

