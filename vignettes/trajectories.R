## ----options, include=F, cache=F, results='hide', message=F----------------

knitr::opts_chunk$set(fig.align="center", cache=FALSE,error=FALSE,
                      fig.width=6,fig.height=6,autodep=TRUE,
                      out.width="600px", out.height="600px",
                      results="markup", echo=TRUE, eval=TRUE)

options(getClass.msg=FALSE)

set.seed(6473) ## for reproducibility


## ---- collapse=F, message=T, eval=F----------------------------------------
#  devtools::install_github("YosefLab/VISION")
#  devtools::install_github("dynverse/dyno")
#  install.packages('tidyverse')

## ---- collapse=F, message=F, warning=F, eval=F-----------------------------
#  library(VISION)
#  library(dyno)
#  library(tidyverse)

## ---- collapse=F, message=T, eval=F----------------------------------------
#  
#  counts = as.matrix(read.table("data/hemato_counts.csv.gz", sep=',', header=T, row.names=1))
#  

## ---- collapse=F, message=T, eval=F----------------------------------------
#  k.genes = read.table("data/bBM.filtered_gene_list.paper.txt.gz", sep=',')[,1]
#  filt.counts = counts[k.genes,]
#  
#  f.genes = VISION:::filterGenesFano(filt.counts)
#  filt.counts = filt.counts[f.genes,]
#  

## ---- collapse=F, message=T, eval=F----------------------------------------
#  
#  scale.factor = median(colSums(counts))
#  expr = apply(counts, 2, function(x) (x * scale.factor) / (sum(x) + 1))
#  filt.expr = log(1 + expr[f.genes,])
#  
#  dataset = wrap_expression(
#    counts = t(filt.counts),
#    expression = t(filt.expr)
#  )
#  
#  model <- infer_trajectory(dataset, "projected_slingshot", verbose=T)
#  
#  model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = dataset$expression)
#  plot_dimred(
#    model,
#    expression_source = dataset$expression
#  )
#  

## ---- collapse=F, eval=F---------------------------------------------------
#  covar = read.table("data/hemato_covariates.txt.gz", sep='\t', header=T, row.names=1)
#  covar = covar[colnames(expr), -1]
#  vis <- Vision(expr,
#              c("data/h.all.v5.2.symbols.gmt"),
#              projection_genes = f.genes,
#              meta=covar,
#              latentTrajectory = model,
#              sig_norm_method="znorm_columns")
#  vis <- addProjection(vis, "MDS", model$dimred[,c(1,2)])
#  
#  vis <- analyze(vis)

## ---- collapse=F, message=T------------------------------------------------
sessionInfo()

