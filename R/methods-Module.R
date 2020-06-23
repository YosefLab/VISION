calcHotspotModules <- function(object, model="danb", tree=F) {
    hotspot <- import("hotspot", convert=F)
    
    if (tree) {
        nwk <- write.tree(object@Tree)
        pyTree <- ete3$Tree(nwk, format = 8L)
        hs <- hotspot$Hotspot(as.data.frame(object@exprData), tree=pyTree, model=model)
    } else {
        hs <- hotspot$Hotspot(as.data.frame(object@exprData), latent=as.data.frame(object@LatentSpace), model=model)
    }
    
    
    hs$create_knn_graph(F, n_neighbors = 30L)
    hs_results <- hs$compute_autocorrelations()
    hs_genes <- hs_results$loc[hs_results$FDR$le(0.05)]$sort_values('Z', ascending=F)$head(1000L)$index
    hs$compute_local_correlations(hs_genes)
    hs$create_modules()
    hs$calculate_module_scores()
    hs_modules <- py_to_r(hs$modules)
    
    
    modules <- object@Modules
    
    for (i in unique(c(hs_modules))) {
      if (i !=-1) {
        names <- dimnames(hs_modules[hs_modules == i])[[1]]
        v <- rep(1, length(names))
        names(v) <- names
        if (tree) {
            new_sig <- Signature(sigDict = v, name=paste("Hotspot Tree Module", i), source="Hotspot", meta=paste("Hotspot on tree, model:", model))
            modules[[paste("HOTSPOT_TREE", i, sep = "_")]] <- new_sig
        } else {
            new_sig <- Signature(sigDict = v, name=paste("Hotspot Module", i), source="Hotspot", meta=paste("Hotspot, model:", model))
            modules[[paste("HOTSPOT", i, sep = "_")]] <- new_sig
        }
      }
    }
    object@Modules <- modules
    
    return(object)
}


