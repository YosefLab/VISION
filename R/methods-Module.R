calcHotspotModules <- function(object, model="normal", tree=FALSE, number_top_genes=1000,
                    num_umi=NULL, min_gene_threshold=20, n_neighbors=NULL,
                    fdr_threshold=0.05) {

    hotspot <- import("hotspot", convert=F)

    workers <- getOption("mc.cores")
    if (is.null(workers)){
        workers <- 1
    }

    exprData = matLog2(object@exprData)

    gene_subset = object@params$latentSpace$projectionGenes

    if (is.null(gene_subset)) {
      gene_subset <- applyFilters(exprData,
            object@params$latentSpace$projectionGenesMethod,
            object@params$latentSpace$threshold, 2)
    }

    exprData = as.data.frame(as.matrix(object@exprData)[gene_subset,])

    # remove genes that do not have any standard deviation
    sds = apply(exprData, 1, sd)
    exprData = exprData[which(sds > 0), ]

    # TODO add UMI support
    if (tree) {
        message("Using Tree")
        ete3 <- import("ete3", convert=F)
        nwk <- write.tree(object@tree)
        pyTree <- ete3$Tree(nwk, format = 8L)
        if (is.null(num_umi)) {
            hs <- hotspot$Hotspot(exprData, tree=pyTree, model=model)
        } else {
            py$umi_df <- r_to_py(num_umi)
            py_run_string("umi_counts = umi_df.iloc[:, 0]")
            hs <- hotspot$Hotspot(exprData, tree=pyTree, model=model, umi_counts=py$umi_counts)
        }
        
    } else {
        if (is.null(num_umi)) {
            hs <- hotspot$Hotspot(exprData, latent=as.data.frame(object@LatentSpace), model=model)
        } else {
            py$umi_df <- r_to_py(num_umi)
            py_run_string("umi_counts = umi_df.iloc[:, 0]")
            hs <- hotspot$Hotspot(exprData, latent=as.data.frame(object@LatentSpace), model=model, umi_counts=py$umi_counts)
        }
    }
    
    if (is.null(n_neighbors)) {
        hs$create_knn_graph(F, n_neighbors = as.integer(object@params$numNeighbors))
    } else {
        hs$create_knn_graph(F, n_neighbors = as.integer(n_neighbors))
    }
  
    hs_results <- hs$compute_autocorrelations(jobs=as.integer(workers))
    hs_genes <- hs_results$loc[hs_results$FDR$le(fdr_threshold)]$sort_values('Z', ascending=F)$head(as.integer(number_top_genes))$index
    hs$compute_local_correlations(hs_genes, jobs=as.integer(workers))
    hs$create_modules(min_gene_threshold=as.integer(min_gene_threshold))
    hs_module_scores <- py_to_r(hs$calculate_module_scores())
    hs_modules <- py_to_r(hs$modules)
    
    
    modules <- list()
    
    for (i in unique(c(hs_modules))) {
      if (i !=-1) {
        names <- dimnames(hs_modules[hs_modules == i])[[1]]
        v <- rep(1, length(names))
        names(v) <- names
        if (tree) {
            new_sig <- Signature(sigDict = v, name=paste("HOTSPOT_TREE", i, sep = "_"), source="Hotspot", meta=paste("Hotspot on tree, model:", model))
            modules[[paste("HOTSPOT_TREE", i, sep = "_")]] <- new_sig
        } else {
            new_sig <- Signature(sigDict = v, name=paste("HOTSPOT", i, sep = "_"), source="Hotspot", meta=paste("Hotspot, model:", model))
            modules[[paste("HOTSPOT", i, sep = "_")]] <- new_sig
        }
      }
    }
    
    colnames(hs_module_scores) <- sort(names(modules))
    object@modData <- modules
    
    if (length(object@sigData) > 0) {
      # Have signatures, let's compute the enrichment
      message("Computing Module-Signature Enrichment")
      object@ModuleSignatureEnrichment <- calc_mod_sig_enrichment(object)
    }
    
    # calculate overlap signatures
    object <- generateOverlapSignatures(object)
    
    # normal analysis
    object <- calcModuleScores(object)
    object <- clusterModScores(object)
    object@Hotspot <- list(hs)
    object@ModuleHotspotScores <- hs_module_scores
    object <- analyzeLocalCorrelationsModules(object, tree)
    
    
    
    return(object)
}




#' calculate module scores
#'
#' For each module-cell pair, compute a score that captures the level of
#' correspondence between the cell and the module.
#'
#' @param object the VISION object
#' @param mod_norm_method Method to apply to normalize the expression matrix
#' before calculating signature scores. Valid options are:
#' "znorm_columns" (default), "none", "znorm_rows", "znorm_rows_then_columns",
#' or "rank_norm_columns"
#' @param mod_gene_importance whether or not to rank each gene's contribution to
#' the overall signature score.  Default = TRUE.  This is used for inspecting
#' genes in a signature in the output report
#' @return the VISION object, with the @ModScores and @ModGeneImportance slots populated
#' @export
calcModuleScores <- function(
    object, mod_norm_method = NULL, mod_gene_importance = TRUE) {
    
    message("Evaluating module scores on cells...\n")
    
    ## override object parameters
    if (!is.null(mod_norm_method)) object@params$modules$modNormMethod <- mod_norm_method
    
    if(is.null(mod_norm_method) && is.null(object@params$modules$modNormMethod)) object@params$modules$modNormMethod <- mod_norm_method <- object@params$signatures$sigNormMethod
    
    if (length(object@modData) == 0) {
      modScores <- matrix(nrow = ncol(object@exprData), ncol = 0,
                          dimnames = list(colnames(object@exprData), NULL)
      )
      object@ModScores <- modScores
      object@ModGeneImportance <- list()
      return(object)
    }
    
    normExpr <- getNormalizedCopySparse(
      object@exprData,
      object@params$modules$modNormMethod
    )
    
    
    modScores <- batchSigEvalNorm(object@modData, normExpr)
    
    if (mod_gene_importance) {
      
      if (is(object@exprData, "sparseMatrix")) {
        modGeneImportance <- evalSigGeneImportanceSparse(
          modScores, object@modData, normExpr
        )
      } else {
        normExprDense <- getNormalizedCopy(
          object@exprData,
          object@params$modules$modNormMethod
        )
        modGeneImportance <- evalSigGeneImportance(
          modScores, object@modData, normExprDense
        )
      }
      
    } else {
      modGeneImportance <- list()
    }
    
    object@ModScores <- modScores
    object@ModGeneImportance <- modGeneImportance
    
    return(object)
}




#' Compute local correlations for all modules
#'
#' This is the main analysis function. For each filtered dataset, a set of
#' different projection onto low-dimensional space are computed, and the
#' consistency of the resulting space with the signature scores is computed
#' to find signals that are captured successfully by the projections.
#' @param object the VISION object
#' @return the VISION object with values set for the analysis results
#' @export
analyzeLocalCorrelationsModules <- function(object, tree=FALSE) {
  
  signatureBackground <- generatePermutationNull(
    object@exprData, object@modData, num = 3000
  )
  
  normExpr <- getNormalizedCopySparse(
    object@exprData,
    object@params$modules$modNormMethod)
  
  if (!tree) {
    message("Computing KNN Cell Graph in the Latent Space...\n")
    weights <- computeKNNWeights(object@LatentSpace, object@params$numNeighbors)
  } else {
    message("Using Tree to compute neighbors...\n")
    weights <- computeKNNWeights(object@tree, object@params$numNeighbors)
  }
  
  message("Evaluating local consistency of modules in latent space...\n")
  
  modConsistencyScores <- sigConsistencyScores(
    weights,
    object@ModScores,
    object@metaData,
    signatureBackground,
    normExpr)
  
  
  modConsistencyScores <- modConsistencyScores[
    colnames(object@ModScores), , drop = FALSE
  ]
  
  
  LocalAutocorrelation <- object@LocalAutocorrelation
  LocalAutocorrelation$Modules = modConsistencyScores
  
  object@LocalAutocorrelation <- LocalAutocorrelation
  
  return(object)
}

calc_mod_sig_enrichment <- function(object) {
    modules <- object@modData
    original_signatures <- object@sigData
    signatures <- list()
    sig_names <- list()
    for (signature in original_signatures) {
        directional <- all(c(1, -1) %in% signature@sigDict)
        if (directional) {
          up <- names(which(signature@sigDict == 1))
          down <- names(which(signature@sigDict == -1))
          
          up_name <- paste(signature@name, "_UP", sep = "")
          down_name <- paste(signature@name, "_DOWN", sep = "")
          
          signatures <- c(signatures, list(up))
          signatures <- c(signatures, list(down))
          
          sig_names <- append(sig_names, up_name)
          sig_names <- append(sig_names, down_name)
        } else {
            signatures <-  c(signatures, list(names(signature@sigDict)))
            sig_names <- append(sig_names, signature@name)
        }
    }
    
    genes <- rownames(object@exprData)
    
    stats <- c()
    p_values <- c()
    for (signature in signatures) {
        set1 <- signature
        stat <- c()
        pval <- c()
        for (module in modules) {
            set2 <- names(module@sigDict)
            results <- calc_set_enrichment(set1, set2, genes)
            stat <- append(stat, results[1])
            pval <- append(pval, results[2])
        }
        stats <- rbind(stats, stat)
        p_values <- rbind(p_values, pval)
    }
    
    mod_names <- c()
    for (module in modules) {
        mod_names <- append(mod_names, module@name)
    }
    
    colnames(stats) <- mod_names
    rownames(stats) <- sig_names
    
    colnames(p_values) <- mod_names
    rownames(p_values) <- sig_names
    
    assignments <- group_modules_enrichment(stats, p_values)
    
    return(list("statistics"=stats, "p_vals"=p_values, "cl"=assignments))
}



calc_set_enrichment <- function(set1, set2, genes) {
    N <- length(genes)
    m <- max(length(set1), length(set2))
    n <- N - m
    k <- min(length(set1), length(set2))
    
    
    o_overlap <- length(intersect(set1, set2))
    e_overlap <- k * (m / N)
    stat <- log(o_overlap / e_overlap)
    if (o_overlap == 0) {
        stat <- 0
    }
    
    p_value <- 1 - phyper(q=o_overlap-1, m, n, k)
    return(c(stat, p_value))
}


group_modules_enrichment <- function(stats, pvals) {
  sigs <- rownames(stats)
  mods <- colnames(stats)
  
  assignments <- as.list(max.col(stats))
  num_cl <- length(unique(assignments <- assignments))
  names(assignments) <- rownames(stats)
  return(assignments)
}


generateOverlapSignatures <- function(object) {
    message("Generating Module Signature Overlaps...\n")
    sigs <- rownames(object@ModuleSignatureEnrichment$statistics)
    mods <- names(object@modData)
    overlap_sigs <- list()
    for (mod in mods) {
      genes <- names(object@modData[[mod]]@sigDict)
      for (sig in sigs) {
        stat <- object@ModuleSignatureEnrichment$statistics[sig, mod]
        pval <- object@ModuleSignatureEnrichment$p_vals[sig, mod]
        if (pval < 0.05) {
          # create overlap signature
          sig_overlap <- object@sigData[[sig]]
          sig_genes <- names(sig_overlap@sigDict)
          overlap_genes <- intersect(genes, sig_genes)
          sig_overlap@sigDict <- sig_overlap@sigDict[overlap_genes]
          new_name <-  paste(sig_overlap@name, "_OVERLAP_", mod, sep="")
          sig_overlap@name <- new_name
          object@modData[[new_name]] <- sig_overlap
        }
      }
    }
    
    return(object)
}




#' Compute Ranksums Test, for all factor meta data.  One level vs all others
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom matrixStats colRanks
#' @importFrom stats setNames
#' @param object the VISION object
#' @param variables which columns of the meta-data to use for comparisons
#' @return the VISION object with the @ClusterComparisons modules slot populated
#' @export
clusterModScores <- function(object, variables = "All") {
  
  message("Computing differential modules and overlaps tests...\n")
  
  modScores <- object@ModScores
  metaData <- object@metaData
  
  metaData <- metaData[rownames(modScores), , drop = FALSE]
  
  if (variables == "All") {
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
  } else {
    if (!all(variables %in% colnames(metaData))) {
      stop("Supplied variable names must be column names of object@metaData")
    }
    clusterMeta <- setNames(variables, variables)
  }
  
  # Comparisons for Modules
  if (ncol(modScores) > 0){
    modScoreRanks <- colRanks(modScores,
                              preserveShape = TRUE,
                              ties.method = "average")
    dimnames(modScoreRanks) <- dimnames(modScores)
  } else {
    modScoreRanks <- modScores
  }
  
  out <- pbmclapply(clusterMeta, function(variable){
    values <- metaData[[variable]]
    var_levels <- levels(values)
    
    result <- lapply(var_levels, function(var_level){
      cluster_ii <- which(values == var_level)
      
      r1 <- matrix_wilcox(modScoreRanks, cluster_ii,
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
  
  object@ClusterComparisons[["Modules"]] <- out
  
  return(object)
  
}


