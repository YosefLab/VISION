simBranchingDataSet <- function(
  scaNCellsPerEdge,
  scaNGroups,
  scaNConst,
  scaNLin,
  scaMumax=100,
  scaSDMuAmplitude=1,
  scaNSigsPerGroup=5,
  vecNormConstExternal=NULL,
  vecDispExternal=NULL,
  vecGeneWiseDropoutRates=NULL,
  matDropoutModelExternal=NULL){

  ####
  # Internal functions
  # Evalute impulse model at time points
  evalImpulse <- function(t,beta,t1,t2,h0,h1,h2){
    return(1/h1* (h0+(h1-h0)*1/(1+exp(-beta*(t-t1))))*
             (h2+(h1-h2)*1/(1+exp(beta*(t-t2)))))
  }

  SCA_MIN_MU <- 10^(-5)
  scaNCells = scaNCellsPerEdge*(scaNGroups-1)

  ### Sample Tree for NGroups
  adjMat <-  matrix(rep(0, scaNGroups^2), scaNGroups, scaNGroups)

  # for simple trajectories
  if (scaNGroups <= 2) {
    adjMat = matrix(c(0,0,1,0), 2,2)
  }
  # Construct higher-order tree
  else {
    children <- c(2:scaNGroups)
    parents <- c(1)
    for (i in 1:(scaNGroups-2)) {
      child <- sample(children,1)
      children <- setdiff(children, child)
      parent <- sample(parents,1)
      adjMat[parent,child] <- 1
      parents <- c(parents, child)
    }
    child = children
    parent <- sample(parents,1)
    adjMat[parent,child] <- 1
  }





  ## assign DE genes to each edge
  arrMuLinIdx <- array(rep(0, scaNGroups^2*scaNLin), c(scaNGroups, scaNGroups, scaNLin))
  scaNcellsPerEdge <- as.integer()
  genesPerEdge <- as.integer(scaNLin/scaNGroups) + 1
  geneSet = c(1:scaNLin)
  edgeMap = c()
  j=1
  for (i in 1:scaNGroups) {
    neighbors<-c(which(adjMat[i,]==1))

    for (neighbor in neighbors) {
      if (j == scaNGroups -1) {
        arrMuLinIdx[i, neighbor, ] <- 1:scaNLin %in% geneSet
        edgeMap = rbind(edgeMap, cbind(rep(i, length(geneSet)), rep(neighbor, length(geneSet))))
      }
      else {
        arrMuLinIdx[i, neighbor,] <- 1:scaNLin %in% geneSet[1:genesPerEdge]
        geneSet <- geneSet[-c(1:genesPerEdge)]
        edgeMap = rbind(edgeMap, cbind(rep(i, genesPerEdge), rep(neighbor, genesPerEdge)))
      }

      j <- j + 1
    }
  }

  dfAnnot <- data.frame(
    cell = c(),
    continuous = c(),
    row.names = c(),
    stringsAsFactors = FALSE
  )

  ## find paths through tree
  findPaths = function(adjMat,
                       paths,
                       curVert,
                       curPath,
                       vecMuLinStart,
                       matMuLinPrev,
                       vecPTPrev,
                       dfAnnot){
    ## Find Neighbors
    neighbors<-c(which(adjMat[curVert,]==1))
    toVisit<-neighbors
    curPath <- c(curPath, curVert)




    #check if end node
    if (length(toVisit) == 0) {
      #end of branch
      i<-length(paths)+1
      paths[paste('PATH', i)] <- list(curPath)
      return(list(paths=paths,
                  matMuLin=matMuLinPrev,
                  dfAnnot=data.frame(cell=colnames(matMuLinPrev),
                                     continuous=vecPTPrev,
                                     row.names = colnames(matMuLinPrev),
                                     stringsAsFactors = FALSE)))
    }

    matMuLinAccum = c()
    dfAnnotAccum = dfAnnot

    while(length(toVisit) > 1) {
      newPath <- curPath
      nextVert <- toVisit[1]

      ## Calculate end gene states
      vecMuLinStart[vecMuLinStart < SCA_MIN_MU] <- SCA_MIN_MU
      vecMuLinEnd <- vecMuLinStart
      idx = which(arrMuLinIdx[curVert, nextVert,]==1)
      vecMuLinEnd[idx] <- vecMuLinStart[idx] * abs(1+rnorm(n=length(idx),
                                                           mean=1,sd=scaSDMuAmplitude))
      vecMuLinEnd[vecMuLinEnd < SCA_MIN_MU] <- SCA_MIN_MU


      vecPT <- seq(0, 1, by=1/(scaNCellsPerEdge-1))
      pathName = paste(tail(c(curPath, nextVert),2), collapse = "")
      cellNames = paste(paste0("CELL_", pathName),seq(1,scaNCellsPerEdge), sep="_")

      matMuLinEdge <- do.call(rbind, lapply(seq_len(scaNLin), function(i) {
        vecMuLinStart[i] + vecPT/max(vecPT)*(vecMuLinEnd[i] - vecMuLinStart[i])
      }))

      rownames(matMuLinEdge) <- vecLinIDs
      colnames(matMuLinEdge) <- cellNames

      vecPT <- (length(curPath)-1)*scaNCellsPerEdge + seq(1,scaNCellsPerEdge)

      res <- findPaths(adjMat,
                       paths,
                       nextVert,
                       newPath,
                       vecMuLinEnd,
                       matMuLinEdge,
                       vecPT,
                       dfAnnot)


      paths = res$paths
      matMuLinAccum = cbind(matMuLinAccum, res$matMuLin)

      if (length(dfAnnotAccum) == 0) {
        dfAnnotAccum = res$dfAnnot
      } else {
        dfAnnotAccum = rbind(dfAnnotAccum, res$dfAnnot)
      }



      toVisit <- toVisit[-1]
    }

    nextVert <- toVisit[1]

    ## Calculate end gene states
    vecMuLinStart[vecMuLinStart < SCA_MIN_MU] <- SCA_MIN_MU
    vecMuLinEnd <- vecMuLinStart
    idx = which(arrMuLinIdx[curVert, nextVert,]==1)
    vecMuLinEnd[idx] <- vecMuLinStart[idx] * abs(1+rnorm(n=length(idx),
                                                         mean=1,sd=scaSDMuAmplitude))
    vecMuLinEnd[vecMuLinEnd < SCA_MIN_MU] <- SCA_MIN_MU

    vecPT <- seq(0, 1, by=1/(scaNCellsPerEdge-1))
    pathName = paste(tail(c(curPath, nextVert),2), collapse = "")
    cellNames = paste(paste0("CELL_", pathName),seq(1,scaNCellsPerEdge), sep="_")

    matMuLinEdge <- do.call(rbind, lapply(seq_len(scaNLin), function(i) {
      vecMuLinStart[i] + vecPT/max(vecPT)*(vecMuLinEnd[i] - vecMuLinStart[i])
    }))

    rownames(matMuLinEdge) <- vecLinIDs
    colnames(matMuLinEdge) <- cellNames
    vecPT <- (length(curPath)-1)*scaNCellsPerEdge + seq(1,scaNCellsPerEdge)


    res <- findPaths(adjMat,
                      paths,
                      nextVert,
                      curPath,
                      vecMuLinEnd,
                      matMuLinEdge,
                      vecPT,
                      dfAnnot)

    matMuLinAccum = cbind(matMuLinAccum,
                          res$matMuLin,
                          matMuLinPrev)

    df =data.frame(cell=colnames(matMuLinPrev),
               continuous=vecPTPrev,
               row.names = colnames(matMuLinPrev),
               stringsAsFactors = FALSE)

    if (length(dfAnnotAccum) == 0) {
      dfAnnotAccum = rbind(res$dfAnnot, df)
    } else {
      dfAnnotAccum = rbind(dfAnnotAccum, res$dfAnnot, df)
    }


    return(list(paths=res$paths,
                matMuLin=matMuLinAccum,
                dfAnnot=dfAnnotAccum))
  }



  ## Init beginning state for Linear DE genes
  vecLinIDs <- paste(paste0(rep("GENE_LIN_",scaNLin),
                            apply(edgeMap, 1, paste, collapse="")),
                      1:scaNLin, sep = "_")


  ## start state for linear DE genes
  vecMuLinStart <- runif(scaNLin)*scaMumax
  vecMuLinStart[vecMuLinStart < SCA_MIN_MU] <- SCA_MIN_MU


  ##
  res <- findPaths(adjMat,
                    list(),
                    1,
                    c(),
                    vecMuLinStart,
                    c(),
                    c(),
                    data.frame())


  ## Draw means from uniform (first half of genes): one mean per gene
  vecConstIDs <- paste0(rep("GENE_CONST_", scaNConst),
                        seq_len(scaNConst))
  vecMuConstHidden <- runif(scaNConst)*scaMumax
  vecMuConstHidden[vecMuConstHidden < SCA_MIN_MU] <- SCA_MIN_MU
  matMuConstHidden <- matrix(vecMuConstHidden,
                             nrow=scaNConst,
                             ncol=scaNCells,
                             byrow=FALSE )

  rownames(matMuConstHidden) <- vecConstIDs
  colnames(matMuConstHidden) <- colnames(res$matMuLin)

  ## Set means for sigs at 1/2 scaMuMax
  vecSigIDs <- paste0(rep("SIG_", scaNGroups*scaNSigsPerGroup),
                      seq_len(scaNGroups*scaNSigsPerGroup))
  vecMuSigHidden <- rep(scaMumax/2, scaNGroups*scaNSigsPerGroup)
  groups <- sapply(strsplit(colnames(res$matMuLin), '_'),
                   function(n) {
                     edges = as.integer(unlist(strsplit(n[2], "")))
                     if (as.integer(n[3]) >= (scaNCellsPerEdge/2)) {
                       return(edges[2])
                     } else{
                       return(edges[1])
                     }
                   })
  matMuSigs = matrix(rep(0, scaNGroups*scaNCells), scaNGroups*scaNSigsPerGroup, scaNCells)


  res$dfAnnot$groups <- groups
  sigs = c()
  for (group in 1:scaNGroups) {
    idx = which(groups == group)
    matMuSigs[((group-1)*scaNSigsPerGroup + 1):((group-1)*scaNSigsPerGroup + scaNSigsPerGroup), idx] <- scaMumax/2
    prefix = paste0('GENE_SIG_', group)
    sigs = c(sigs,
             paste(prefix, 1:scaNSigsPerGroup, sep='_'))
  }
  dfSigs <- data.frame(groups=1:scaNGroups, t(matrix(sigs, scaNSigsPerGroup, scaNGroups)))
  rownames(matMuSigs) = sigs
  matMuHidden <- do.call(rbind, list(matMuConstHidden,
                                     res$matMuLin,
                                     matMuSigs))



  ## Draw dispersion

  vecDispHidden <- 3 + rnorm(n = dim(matMuHidden)[1], mean = 0, sd = 0.5)
  vecDispHidden[vecDispHidden < 0.05] <- 0.05

  matDispHidden <- matrix(vecDispHidden,nrow=dim(matMuHidden)[1],
                          ncol=dim(matMuHidden)[2], byrow=FALSE)
  rownames(matDispHidden) <- rownames(matMuHidden)
  colnames(matDispHidden) <- colnames(matMuHidden)


  ## add noise - draw from negative binomial
  message("Simulate negative binomial noise")
  matSampledDataHidden <- do.call(rbind, lapply(
    seq_len(nrow(matMuHidden)), function(gene){
      sapply(seq_len(scaNCells), function(cell){
        rnbinom(n=1, mu=matMuHidden[gene,cell], size=vecDispHidden[gene])
      })
    }))
  rownames(matSampledDataHidden) <- rownames(matMuHidden)
  colnames(matSampledDataHidden) <- colnames(matMuHidden)


  ### 3. Apply drop out
  message("Simulate drop-out")
  ## Generate underlying drop-out rate matrix
  if(!is.null(matDropoutModelExternal) & !is.null(vecGeneWiseDropoutRates)) {
    stop("Supply either matDropoutModelExternal",
         " or vecGeneWiseDropoutRates.")
  }
  if(!is.null(matDropoutModelExternal)){
    lsDropoutRatesHidden <- lapply(seq_len(scaNCells), function(cell){
      decompressDropoutRateByCell(
        vecDropModel=matDropoutModelExternal[cell,],
        vecMu=matMuHidden[,cell],
        matPiConstPredictors=NULL,
        lsDropModelGlobal=list(strDropModel = "logistic_ofMu",
                               scaNumCells = scaNCells))
    })
    matDropoutRatesHidden <- do.call(cbind, lsDropoutRatesHidden)
    rownames(matDropoutRatesHidden) <- rownames(matMuHidden)
    colnames(matDropoutRatesHidden) <- colnames(matMuHidden)
  } else if(!is.null(vecGeneWiseDropoutRates)) {
    matDropoutRatesHidden <- matrix(
      vecGeneWiseDropoutRates,
      nrow = length(vecGeneWiseDropoutRates),
      ncol = dim(matMuHidden)[2], byrow = FALSE)
  } else {
    stop("Supply either matDropoutModelExternal",
         " or vecGeneWiseDropoutRates.")
  }
  ## Draw drop-outs from rates: Bernoulli experiments
  lsDropoutsHidden <- lapply(
    seq_len(nrow(matDropoutRatesHidden)), function(i){
      rbinom(n=rep(1, dim(matDropoutRatesHidden)[2]),
             size=rep(1, dim(matDropoutRatesHidden)[2]),
             prob=matDropoutRatesHidden[i,])
    })
  matDropoutsHidden <- do.call(rbind, lsDropoutsHidden)
  rownames(matDropoutsHidden) <- rownames(matMuHidden)
  colnames(matDropoutsHidden) <- colnames(matMuHidden)

  ## Create observed data: merge hidden data with drop-outs
  matSampledDataObserved <- matSampledDataHidden
  matSampledDataObserved[matDropoutsHidden==1] <- 0

  ## Round to counts
  matSampledCountsObserved <- round(matSampledDataObserved)


  roots = res$dfAnnot$cell[which(res$dfAnnot$continuous == 1 & res$dfAnnot$groups == 1)]

  scaNGroups = max(dfSigs$groups)
  groups = 1:scaNGroups
  nGenes = nrow(matMuHidden)
  startGene = nGenes - scaNGroups*(dim(dfSigs)[2]-1)

  names = rbind(unlist(sapply(1:scaNGroups, function(i) {
    rep(i, dim(dfSigs)[2] -1)
  },
  simplify = FALSE)))


  names = paste0('GROUP_', names)
  sign = rep('Plus', length(names))

  sigs = cbind(names,
                sign,
                sigs)


  return(list(counts=matSampledCountsObserved,
              dfAnnot=res$dfAnnot,
              adjMat=adjMat,
              arrSignature=arrMuLinIdx,
              matMu=matMuHidden,
              matDisp=matDispHidden,
              roots=roots,
              dfSigs=dfSigs,
              housekeeping = vecConstIDs,
              sigs = sigs))
}


genRandomSigs <- function(expr, nGroups, nGenesPerSig) {
  names = rbind(unlist(sapply(1:nGroups, function(i) {
    rep(i, nGenesPerSig)
  },
  simplify = FALSE)))
  names = paste0('GROUP_', names)
  sign = rep('Plus', length(names))
  sigs = sample(rownames(expr), size=length(names), replace=TRUE)
  
  return(cbind(names,
               sign,
               sigs))
}
