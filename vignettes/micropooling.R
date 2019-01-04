## ----options, include=F, cache=F, results='hide', message=F----------------

knitr::opts_chunk$set(fig.align="center", cache=FALSE,error=FALSE,
                      fig.width=6,fig.height=6,autodep=TRUE,
                      out.width="600px", out.height="600px",
                      results="markup", echo=TRUE, eval=TRUE)

options(getClass.msg=FALSE)

set.seed(6473) ## for reproducibility


## ---- collapse=F, message=T, eval=F----------------------------------------
#  devtools::install_github("YosefLab/VISION")

## ---- collapse=F, message=F, warning=F, eval=F-----------------------------
#  library(VISION)
#  
#  counts = as.matrix(read.table("data/hemato_counts.csv.gz", sep=',', header=T, row.names=1))
#  
#  # compute scaled counts
#  scale.factor = median(colSums(counts))
#  scaled.counts = t(t(counts) / colSums(counts)) * scale.factor
#  
#  # perform preliminary Fano filtering to determing projection genes, as usual
#  f.genes = VISION:::filterGenesFano(scaled.counts)
#  
#  # read in meta data
#  meta = read.table("data/hemato_covariates.txt.gz", sep='\t', header=T, row.names=1)
#  meta = meta[colnames(scaled.counts), -1]
#  
#  vis <- Vision(scaled.counts,
#                c("data/h.all.v5.2.symbols.gmt"),
#                pool=T,
#                cellsPerPartition=5,
#                projection_genes = f.genes,
#                meta=meta)
#  
#  vis <- analyze(vis)
#  
#  viewResults(vis)

## ---- collapse=F, message=F, warning=F, eval=F-----------------------------
#  
#  cellsPerPartition = 5
#  
#  meta$GRvsER = as.factor(sapply(meta$ct, function(x) ifelse((x == "BA" || x == "ER" || x == "MK"), "ErLin", "GrLin")))
#  
#  meta.var = meta[,"GRvsER", drop=F]
#  all.pools = sapply(levels(meta.var$GRvsER), function(x) {
#      print(x)
#      cells = rownames(meta.var)[meta.var[,1] == x]
#      pools = VISION:::applyMicroClustering(scaled.counts[, cells],
#                                cellsPerPartition=cellsPerPartition,
#                                filterInput = f.genes,
#                                filterThreshold = 0.1,
#                                preserve_clusters = NULL,
#                                latentSpace = matrix(NA, 1, 1))
#  
#      nn <- sapply(1:length(pools), function(y) paste0(x, ".microcluster", y))
#  
#      names(pools) <- nn
#      return(pools)
#  
#    })
#  
#  all.pools = unlist(all.pools, recursive=F)
#  
#  vis <- Vision(scaled.counts,
#                c("data/h.all.v5.2.symbols.gmt"),
#                projection_genes = f.genes,
#                meta=meta,
#                pools = all.pools)
#  
#  vis <- analyze(vis)
#  

## ---- collapse=F, message=T------------------------------------------------
sessionInfo()

