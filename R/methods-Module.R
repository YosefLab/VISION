calcHotspotModules <- function(object, model="danb", tree=F, genes=1000, num_umi=NULL, min_gene_threshold=20, n_neighbors=NULL, jobs=1) {
    hotspot <- import("hotspot", convert=F)
    # TODO add UMI support
    if (tree) {
        message("Using Tree to compute knn")
        ete3 <- import("ete3", convert=F)
        nwk <- write.tree(object@Tree)
        pyTree <- ete3$Tree(nwk, format = 8L)
        if (is.null(num_umi)) {
            hs <- hotspot$Hotspot(as.data.frame(object@exprData), tree=pyTree, model=model)
        } else {
            hs <- hotspot$Hotspot(as.data.frame(object@exprData), tree=pyTree, model=model, umi_counts=num_umi)
        }
        
    } else {
        if (is.null(num_umi)) {
            hs <- hotspot$Hotspot(as.data.frame(object@exprData), latent=as.data.frame(object@LatentSpace), model=model)
        } else {
            hs <- hotspot$Hotspot(as.data.frame(object@exprData), latent=as.data.frame(object@LatentSpace), model=model, umi_counts=num_umi)
        }
    }
    
    if (is.null(n_neighbors)) {
        hs$create_knn_graph(F, n_neighbors = as.integer(object@params$numNeighbors))
    } else {
      hs$create_knn_graph(F, n_neighbors = as.integer(n_neighbors))
    }
    
    hs_results <- hs$compute_autocorrelations(jobs=as.integer(jobs))
    hs_genes <- hs_results$loc[hs_results$FDR$le(0.05)]$sort_values('Z', ascending=F)$head(as.integer(genes))$index
    hs$compute_local_correlations(hs_genes, jobs=as.integer(jobs))
    hs$create_modules(min_gene_threshold=as.integer(min_gene_threshold))
    # hs$calculate_module_scores()
    hs_modules <- py_to_r(hs$modules)
    #scores <- as.matrix(py_to_r(hs$module_scores))
    
    modules <- object@modData
    # mod_scores <- object@ModScores
    # 
    # if (all(is.na(mod_scores))) {
    #     mod_scores <- scores
    # } else {
    #     mod_scores <- cbind(mod_scores, scores)
    # }
    # 
    # for (i in unique(c(hs_modules))) {
    #   if (i !=-1) {
    #     if (tree) {
    #       name <- paste("HOTSPOT_TREE", i, sep = "_")
    #     } else {
    #       name <- paste("HOTSPOT", i, sep = "_")
    #     }
    #     colnames(mod_scores) <- append(colnames(mod_scores), name)
    #   }
    # }
    # 

    
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
    
    object@modData <- modules
    object <- calcModuleScores(object)
    
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
    weights <- computeKNNWeights(object@Tree, object@params$numNeighbors)
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


