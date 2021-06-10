#' Perform Hotspot analysis on Vision Object
#'
#' @param object Vision Object
#' @param model model argument for Hotspot, one of \itemize{
#' \item normal
#' \item danb
#' \item bernoulli
#' \item none 
#' }
#' @param tree whether to use tree as latent space. If TRUE, object should have
#' a tree slot.
#' @param number_top_genes Hotspot argument for number of genes to consider
#' @param num_umi optional dataframe containing umi counts in first column for
#'  barcodes
#' @param min_gene_threshold minimum number of genes in Hotspot module
#' @param n_neighbors number of neighbors to consider in latent space
#' @param autocorrelation_fdr threshold for significance for genes autocorr
#' @param clustering_fdr threshold for significance for clustering modules
#' @param logdata boolean, log the expression data, avoid for danb
#' Populates the modData, HotspotModuleScores, ModuleSignatureEnrichment
#' and HotspotObject slots of object, as well as recalculates signature scores
#' for new modules.
#' @return the modified Vision object
#' 
#' @export
runHotspot <- function(object, model="normal", tree=FALSE, 
                               number_top_genes=1000, num_umi=NULL, 
                               min_gene_threshold=20, n_neighbors=NULL,
                               autocorrelation_fdr=0.05, clustering_fdr=0.5, logdata=FALSE) {

    # Init Hotspot
    hs <- hsInit(object, model, tree, num_umi, logdata)
    # Init Hotspot KNN
    hs <- hsCreateKnnGraph(hs, object, n_neighbors=n_neighbors)
    # perform Hotspot analysis and store results in R
    hs_genes <- hsComputeAutoCorrelations(hs, number_top_genes=number_top_genes, autocorrelation_fdr=autocorrelation_fdr)
    # Compute localcorr
    hs <- hsComputeLocalCorrelations(hs, hs_genes)
    # Calculate Hotspot Module Scores for informative genes
    hs <- hsCalculateModuleScores(hs, min_gene_threshold, clustering_fdr)
    # Cluster Hotspot modules and perform Vision based analysis on HS Modules and 
    object <- analyzeHotspotObjectVision(object, hs, tree)
    
    return(object)
}


#' Init Hotspot object from Vision Object
#' 
#' @param object the Vision Object
#' @param model the model for Hotspot (ie "normal", "danb"...)
#' @param tree boolean, whether to use the tree as ls
#' @param num_umi df of barcodes x num_umi
#' @param logdata boolean, log the expression data, avoid for danb
#' 
#' @return the Hotspot object
#' 
#' @export
hsInit <- function(object, model="normal", tree=F, num_umi=NULL, logdata=FALSE) {
  hotspot <- import("hotspot", convert=F)
  
  workers <- getOption("mc.cores")
  if (is.null(workers)){
    workers <- 1
  }
  
  if (!logdata) {
    # Don't take the log
    exprData = object@exprData
  } else {
    # take the log2 otherwise
    exprData = matLog2(object@exprData)
  }
  
  gene_subset <- object@params$latentSpace$projectionGenes
  
  if (any(is.na(gene_subset))) {
    gene_subset <- applyFilters(exprData,
                                object@params$latentSpace$projectionGenesMethod,
                                object@params$latentSpace$threshold, 2)
  }
  
  exprData = as.data.frame(as.matrix(exprData)[gene_subset,])
  
  # remove genes that do not have any standard deviation
  sds = apply(exprData, 1, sd)
  exprData = exprData[which(sds > 0), ]
  
  # generate the Hotspot object in python, potentially using the tree
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
  
  return(hs)
}


#' Init KNN graph in Hotspot object
#' 
#' @return the Hotspot object with KNN initialized
#' 
#' @export
hsCreateKnnGraph <- function(hs, object, n_neighbors=NULL) {
  # create knn graph, specify nn or use object default
  if (is.null(n_neighbors)) {
    hs$create_knn_graph(F, n_neighbors = as.integer(object@params$numNeighbors))
  } else {
    hs$create_knn_graph(F, n_neighbors = as.integer(n_neighbors))
  }
  return(hs)
}


#' Compute Hotspot auto correlations
#' 
#' @param hs the Hotspot object
#' @param number_top_genes Hotspot argument for number of genes to consider
#' @param autocorrelation_fdr threshold for significance for genes autocorr
#' @return list of HS genes
#' 
#' @export
hsComputeAutoCorrelations <- function(hs, number_top_genes=1000, autocorrelation_fdr=0.05) {
  workers <- getOption("mc.cores")
  if (is.null(workers)){
    workers <- 1
  }
  
  hs_results <- hs$compute_autocorrelations(jobs=as.integer(workers))
  hs_genes <- hs_results$loc[hs_results$FDR$le(autocorrelation_fdr)]$sort_values('Z', ascending=F)$head(as.integer(number_top_genes))$index
  return(hs_genes)
}


#' Interface function to compute local correlations for Hotspot
#' Warning: modifies the hs argument
#' @param hs the Hotspot object
#' @param hs_genes Hotspot genes
#' @return the populated hs object
#' 
#' @export
hsComputeLocalCorrelations <- function(hs, hs_genes) {
  workers <- getOption("mc.cores")
  if (is.null(workers)){
    workers <- 1
  }
  hs$compute_local_correlations(hs_genes, jobs=as.integer(workers))
  return(hs)
}


#' Analyze a Hotspot object using built in methods such
#' such as local correlation, signature overlap, etc.
#' Necessary to run this function for Hotspot functionality in viewer to work.
#'
#' @param object the VISION object
#' @param hs the Hotspot python object loaded by Reticulate
#' @param tree whether to use tree as latent space. If TRUE, object should have a tree
#' 
#' @return the modified VISION object with the following slots filled:
#' Populates the modData, HotspotModuleScores, ModuleSignatureEnrichment
#' and HotspotObject slots of object, as well as recalculates signature scores
#' for new modules.
#' 
#' @export
analyzeHotspotObjectVision <- function(object, hs, tree=FALSE) {
    hs_module_scores <- hs$module_scores
    hs_modules <- hs$modules
    
    # handle annoying reticulate conversion issues when writing to a file
    if (!is.data.frame(hs_module_scores)) {
        hs_module_scores <- py_to_r(hs_module_scores)
    }
    
    if (!is.array(hs_modules)) {
        hs_modules <- py_to_r(hs_modules)
    }
    
    modules <- list()
    
    # add the modules with name scheme HOTSPOT_{TREE if using TREE}_#
    model <- hs$model
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
    
    # store module scores
    colnames(hs_module_scores) <- sort(names(modules))
    object@modData <- modules
    
    if (length(object@sigData) > 0) {
      # Have signatures, let's compute the enrichment
      message("Computing Module-Signature Enrichment")
      object@ModuleSignatureEnrichment <- calc_mod_sig_enrichment(object)
    }
    
    # calculate overlap signatures
    object <- generateOverlapSignatures(object)
    
    # re run normal analysis
    object <- calcModuleScores(object)
    object <- clusterModScores(object)
    
    object@ModuleHotspotScores <- hs_module_scores
    object <- analyzeLocalCorrelationsModules(object, tree)
    
    # save the Hotspot object
    object <- addHotspotToVision(object, hs)
    
    return(object)
}


#' Create Hotspot Modules and calculate module scores given a HS object
#' with local correlations already calculated
#' 
#' @param hs the Hotspot object, must have ran compute_local_correlations already
#' @param min_gene_threshold min genes per module
#' @param clustering_fdr p value for clustering genes
#' @return the modified hs object
#' 
#' @export
hsCalculateModuleScores <- function(hs, min_gene_threshold=20, clustering_fdr=0.5, plot=F) {
  hs$create_modules(min_gene_threshold=as.integer(min_gene_threshold), fdr_threshold=clustering_fdr)
  hs$calculate_module_scores()
  
  if (plot) {
    draw_hotspot_heatmap(hs)
  }
  
  return(hs)
}


#' Add HS python obj to vision OBJECT
#' 
#' @param object Vision object
#' @param hs python hs object
#' @return Vision object with hs populated
#' 
#' @export
addHotspotToVision <- function(object, hs) {
    # save the Hotspot object
    pickle <- import("pickle", convert=F)
    py$hs <- hs
    py$pickle <- pickle
    py_run_string("hs_byte_array = bytearray(pickle.dumps(hs))")
    hs_pickled_r <- as.raw(py$hs_byte_array)
    object@Hotspot <- hs_pickled_r
    
    return(object)
}


#' Calculate module scores (signature scores but on the modules)
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
#' 
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
#' 
#' @param object the VISION object
#' @param tree whether to use the tree object as latent space for neighbors
#' @return the VISION object with values set for the analysis results
#' 
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


#' Computes the hypergeometric overlap test for modules and signatures
#' 
#' @param object the Vision object.
#' @param skip_down whether to ignore down signatures in overlap
#' @return list(statistic values, p values, clusters of signatures)
#' 
#' @export
calc_mod_sig_enrichment <- function(object, skip_down=TRUE) {
    modules <- object@modData
    original_signatures <- object@sigData
    signatures <- list()
    sig_names <- list()
    for (signature in original_signatures) {
        # calculate enrichment for both the up signal and down signal signature genes
        directional <- all(c(1, -1) %in% signature@sigDict)
        down_reg <- all(signature@sigDict == -1)
        if (skip_down && down_reg) {
          next
        }
        if (directional) {
          up <- names(which(signature@sigDict == 1))
          down <- names(which(signature@sigDict == -1))
          
          up_name <- paste(signature@name, "_UP", sep = "")
          down_name <- paste(signature@name, "_DOWN", sep = "")
          
          signatures <- c(signatures, list(up))
          sig_names <- append(sig_names, up_name)
          
          if (!skip_down) {
            signatures <- c(signatures, list(down))
            sig_names <- append(sig_names, down_name)
          }
         
        } else {
            signatures <-  c(signatures, list(names(signature@sigDict)))
            sig_names <- append(sig_names, signature@name)
        }
    }
    
    genes <- rownames(object@exprData)
    
    # calculate the enrichments
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
    
    # group the signatures
    assignments <- group_modules_enrichment(stats, p_values)
    
    return(list("statistics"=stats, "p_vals"=p_values, "cl"=assignments))
}


#' Calculate the hypergeometric enrichment for two sets from a population
#' Statisic = log (observed overlap / expected overlap)
#' P value = 1- hypergeometric(observed overlap -1, max(|set1|, |set2|), |genes| - |set1|, min(|set1|, |set2|))
#' 
#' @param set1
#' @param set2
#' @param genes the population
#' @return c(statistic, p value)
#' 
#' @export
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


#' Make the clusters for the modules by enrichment.
#' For now we just assign each signature to each cluster, could filter to only include once,
#' so that each one appears in the modules x sigs table.
#' 
#' @param stats overlap stats from calc_set_enrichment
#' @param pvals overlap p values from calc_set_enrichment
#' @return assignments of each signature to each module
#' 
#' @export
group_modules_enrichment <- function(stats, pvals) {
  sigs <- rownames(stats)
  mods <- colnames(stats)
  
  assignments <- as.list(max.col(stats))
  num_cl <- length(unique(assignments <- assignments))
  names(assignments) <- rownames(stats)
  return(assignments)
}


#' Generates signature objects for the overlap sets between modules and signatures
#' @param object the Vision object
#' @return Vision Object, populates the modData slot with overlap signatures.
#' 
#' @export
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
#' 
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


#' Load in an existing Hotspot object from bytes or a file
#'
#' @param file optional path to an existing file containing pickled bytes (format 0)
#' @param bytes optional R character vector of bytes created by calcHotspotModules
#' @return an externalptr to a Hotspot Object loaded in the R reticulate session
#' 
#' @export
loadHotspotObject <- function(file=NULL, bytes=NULL) {
  hotspot <- import("hotspot", convert=F)
  pickle <- import("pickle", convert=F)
  py$pickle <- pickle
  
  if (!is.null(file)) {
    # load file
    hs <- py_load_object(file)
    return(hs)
  } else if (!is.null(bytes)) {
    py$hs_bytes <- r_to_py(bytes)
    py_run_string("hs = pickle.loads(hs_bytes)")
    return(py$hs)
  } else {
    return(NULL)
  }
}


#' Save bytes in the Hotspot object slot to a file
#' @param path the file path
#' @param bytes the raw bytes, like in the Hotspot slot of a VISION Object
#' 
#' @export
saveHSBytestToPickle <- function(path, bytes) {
  py_save_object(obj=loadHotspotObject(bytes=bytes), path)
}


#' Add custom tree based neighbor and weights to a Hotspot object
#' 
#'  @param tree object of class phylo
#'  @param the Hotspot object to add the nw to
#'  @param minSize the minimum number of neighbors of the node
#'  @return the Hotspot object
#'  
#'  @export
lcaBasedHotspotNeighbors <- function(tree, hotspot, minSize=20) {
  tips <- tree$tip.label
  nTips <- length(tips)
  neighbors <- data.frame(t(matrix(seq_len(nTips) -1, ncol = nTips, nrow= nTips)))
  rownames(neighbors) <- tips
  
  
  weights <- data.frame(matrix(0, ncol = nTips, nrow= nTips))
  for (tip in seq_len(nTips)) {
    my_neighbors <- minSizeCladeNeighbors(tree, tip, minSize)
    
    weights[tip, my_neighbors] <- 1
  }
  
  neighbors_no_diag <- data.frame(matrix(ncol = nTips -1, nrow= nTips))
  weights_no_diag <- data.frame(matrix(ncol = nTips -1, nrow= nTips))
  
  for (tip in seq_len(nTips)) {
    neighbors_no_diag[tip, ] <- neighbors[tip, -tip]
    weights_no_diag[tip, ] <- weights[tip, -tip]
  }
  
  rownames(neighbors_no_diag) <- tips
  rownames(weights_no_diag) <- tips
  
  colnames(neighbors_no_diag) <- seq_len(nTips-1) - 1
  colnames(weights_no_diag) <- seq_len(nTips-1) - 1
  return(list("neighbors"=neighbors_no_diag, "weights"=weights_no_diag))
}


#' Draw Modules Heatmap (Gene x Gene)
#' 
#' @param hs the Hotspot Object
#' @param palette palette
#' @export
draw_hotspot_heatmap <- function(hs, palette = paletteer_d("ggsci::default_nejm")) {
  scipy_hierarchy <- import("scipy.cluster.hierarchy", convert=F)
  np <- import("numpy")
  linkage <- hs$linkage
  py_dend = scipy_hierarchy$dendrogram(linkage)
  lcz = hs$local_correlation_z
  gene_order = colnames(lcz)[np$array(py_dend$leaves)+1]
  col_mapping = c()
  module_to_col = list()
  unique_mods = as.character(unique(hs$modules))
  example_genes = list()
  col_mapping[["-1"]] = "#ffffff"
  for (i in 1:length(unique_mods)) {
    mod = unique_mods[[i]]
    if (mod != -1) {
      col_mapping[[mod]] = palette[[i]] 
      module_to_col[[mod]] = palette[[i]]
      example_genes[[mod]] = sample(hs$modules[hs$modules == mod], 5)
    }
  }
  print(module_to_col)
  print(col_mapping[order(names(col_mapping))])
  print(example_genes)
  modules = data.frame("module" = hs$modules)
  modules$module = modules$module
  ha = rowAnnotation(df = modules,
                     col = col_mapping,
                     simple_anno_size = unit(0.5, "in"))
  ht = Heatmap(as.matrix(lcz), name = "mat",
               show_row_names=F, show_column_names=F, show_row_dend=F, show_column_dend=F,
               row_order = gene_order, column_order=gene_order,
               right_annotation=ha,       
               column_names_gp = gpar(fontsize = c(1)),
               width = unit(8, "in"), height = unit(8, "in"))
  draw(ht)
}
