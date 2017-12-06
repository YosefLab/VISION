#' Computes path statistics for a TreeProjection starting from a given
#' root or search through all vertices
#'
#' @param projection TreeProjection object
#' @param root cell to use for the beginning of the pseudotime projection
#' default: NULL
#' @param k numeric max number of vertices in path. k=0 => full paths
#' default: 0
#' @return list:
#' \itemize{
#'     \item path i: a list containing the cells and pseudotime coordinates
#'     for each path originating from root.
#' }
#' @export
testProjectionLineages <- function(projection, fcounts, root=NULL, k=0) {
  #Check to run from given root
  if (is.null(root)) {
    print("TBD finish unsupervised search")
    return(NULL)
  }


  pt <- pseudotimeProjection(projection, root)

  #test each path
  lsObjLP <- lapply(pt, function(path) {
    objLP <- testPath(path$data, fcounts)

  })
}

testPathLineages <- function(projection, paths, root=NULL, k=0) {
  #Check to run from given root
  if (is.null(root)) {
    print("TBD finish unsupervised search")
    return(NULL)
  }


  #test each path
  lsObjLP <- lapply(paths, function(path) {
    objLP <- testPath(path$data, projection@pData)

  })


  return(lsObjLP)
}
#' Filter Lineages
#' @param lsObjLP list LineagePUlse objects
#' @return list: containg the signifigant genes corresponding to LineagePulse obj
#' \itemize{
#'     \item sigGenes i
#' }
#' @export
testSigGenes <- function(lsObjLP) {
  #filter results
  lsSigGenes <- lapply(lsObjLP, function(objLP) {
    df <- objLP$dfResults
    sigGenes <- row.names(df[which(df$padj <= 0.05),])
    return(sigGenes)

  })
  return(lsSigGenes)
}


#' Computes pseudotime projection for a given TreeProjection and root.
#'
#' @param projection TreeProjection object
#' @param root cell to use for the beginning of the pseudotime projection
#' @return list:
#' \itemize{
#'     \item path i: a list containing the cells and pseudotime coordinates
#'     for each path originating from root.
#' }
#' @export
pseudotimeProjection <- function(projection, root) {
  # find root cell tree vertex
  rootIdx <- match(root, colnames(projection@pData))
  rootPos <- projection@edgePos[rootIdx]
  if (rootPos >= 0.5) {
    rootVert <- projection@edgeAssoc[2, rootIdx]
  }
  else {
    rootVert <- projection@edgeAssoc[1, rootIdx]
  }
  paths <- list()
  notVisited <- c(1:length(projection@adjMat[1,]))



  findPaths <- function(projection,
                       paths,
                       curVert,
                       prevVert,
                       curPath) {

    neighbors <- c(which(projection@adjMat[curVert,]==1))
    toVisit <- neighbors[which(neighbors != prevVert)]
    curPath$verts <- c(curPath$verts, curVert)



    #check if end node
    if (length(toVisit) == 0) {
      #end of branch
      i <- length(paths)+1
      paths[paste('path', i)] <- list(curPath)
      return(paths)
    }



    while (length(toVisit) > 1) {
      # create new branch
      newPath <- curPath

      nextVert <- toVisit[1]
      #find cells on path
      if (nextVert > curVert) {
        edge <- c(curVert, nextVert)
      }
      else {
        edge <- c(nextVert, curVert)
      }

      cellIdx <- which(apply(projection@edgeAssoc, 2, function(i) all(edge == i)))
      cells <- colnames(projection@pData)[cellIdx]

      edgeVec <- projection@vData[,edge[2]] - projection@vData[,edge[1]]

      cellPos <- projection@edgePos[cellIdx]

      relProj <- t(t(matrix(rep(edgeVec, length(cellPos)),
                           nrow=length(edgeVec),
                           ncol=length(cellPos)))*cellPos)


      #cellProj = relProj + projection@vData[,edge[1]]
      cellProj <- projection@spatialPos[,cellIdx]

      #cellDist = curPath$distance + cellPos*sqrt(sum(edgeVec^2))
      if (is.null(dim(cellProj))) {
        cellDist <- curPath$distance + sqrt(sum((cellProj-edgeVec)^2))
      } else {
        cellDist <- curPath$distance + sqrt(colSums((cellProj-edgeVec)^2))
      }

      #add cells to path
      newPath$cells <- c(curPath$cells, cells)
      newPath$cellProj <- cbind(curPath$cellProj, cellProj)
      newPath$cellDist <- c(curPath$cellDist, cellDist)
      newPath$distance <- curPath$distance + sqrt(sum(edgeVec^2))


      paths <- findPaths(projection,
                        paths,
                        nextVert,
                        curVert,
                        newPath)

      toVisit <- toVisit[-1]


    }

    nextVert <- toVisit[1]
    #find cells on path
    if (nextVert > curVert) {
      edge <- c(curVert, nextVert)
    }
    else {
      edge <- c(nextVert, curVert)
    }

    cellIdx <- which(apply(projection@edgeAssoc, 2, function(i) all(edge == i)))
    cells <- colnames(projection@pData)[cellIdx]

    edgeVec <- projection@vData[,edge[2]] - projection@vData[,edge[1]]

    cellPos <- projection@edgePos[cellIdx]

    relProj <- t(t(matrix(rep(edgeVec, length(cellPos)),
                         nrow=length(edgeVec),
                         ncol=length(cellPos)))*cellPos)


    #cellProj = relProj + projection@vData[,edge[1]]
    cellProj <- projection@spatialPos[,cellIdx]

    #cellDist = curPath$distance + cellPos*sqrt(sum(edgeVec^2))
    if (is.null(dim(cellProj))) {
      cellDist <- curPath$distance + sqrt(sum((cellProj-edgeVec)^2))
    } else {
      cellDist <- curPath$distance + sqrt(colSums((cellProj-edgeVec)^2))
    }


    #add cells to path
    curPath$cells <- c(curPath$cells, cells)
    curPath$cellProj <- cbind(curPath$cellProjs, cellProj)
    curPath$cellDist <- c(curPath$cellDist, cellDist)
    curPath$distance <- curPath$distance + sqrt(sum(edgeVec^2))


    paths <- findPaths(projection,
                      paths,
                      nextVert,
                      curVert,
                      curPath)

    return(paths)
  }


  curPath <- list(cells=c(), cellProj=c(), cellDist= c(), distance=0, verts=c())
  paths <- findPaths(projection,
                    paths,
                    rootVert,
                    rootVert,
                    curPath)

  maxDist <- 0
  minDist <- min(paths$`path 1`$cellDist)

  for (path in paths) {
    pathMax <- max(path$cellDist)
    pathMin <- min(path$cellDist)
    if (pathMax > maxDist) {
      maxDist <- pathMax
    }
    if (pathMin < minDist) {
      minDist <- pathMin
    }
  }
  distRange <- maxDist - minDist

  pathReturn <- list()
  for (path in paths) {
    i <- length(pathReturn) + 1
    pts <- (path$cellDist - minDist) / distRange
    df <- data.frame(cell = path$cells, continuous=pts)
    rownames(df) <- df$cell
    pathReturn[paste('path', paste(path$verts, collapse = "_"))] <- list(list(verts=path$verts, data=df))
  }
  return(pathReturn)
}

#' Utility function for computing the closest vertex for each cell in a projection
#'
#' @param projection TreeProjection object
#' @return numeric: closest vertex for each cell in projection
getVerts <- function(projection) {
  assignments <- rep(1, length(projection@edgePos))
  idx <- projection@edgePos >= 0.5
  assignments[idx] <- 2
  verts <- sapply(1:length(assignments), function(i) { projection@edgeAssoc[assignments[i],i]})
  names(verts) <- colnames(projection@pData)
  return(verts)
}


#' Utility function for finding the vertex for a given cell in a projection
#'
#' @param projection TreeProjection object
#' @param cell character name of cell
#' @return numeric: closest vertex for given cell in projection
getCellVert <- function(projection, cell) {
  cellIdx = match(cell, colnames(projection@pData))
  verts = getVerts(projection)
  return(verts[cellIdx])
}


#' Utility function for running LineagePulse with default arguments
#'
#' @param path data.frame containing cells names and continous pseudotime
#' @param counts matrix genes X cells
#' @return numeric: closest vertex for given cell in projection
testPath <- function(path, counts) {
  counts = counts[, levels(path$cell)]
  objLP <- runLineagePulse(
    counts = counts,
    dfAnnotation = path,
    strMuModel = "splines",
    scaDFSplinesMu = 3)
  return(objLP)
}


#' Utility function converting slingshot trajecories to LineagePulse input
#' arguments.
#'
#' @param ptSLing matrix numeric cells X curves
#' @return data.frame: index cells columns curves
ptSling.to.ptLin <- function(ptSling) {
  if (is.null(dim(ptSling))) {
    idx = order(ptSling, na.last=NA)
    return(data.frame(cells=names(ptSling[idx]),
                      continuous=ptSling[idx]))
  }

  return(lapply(1:dim(ptSling)[2], function(i) {
    idx = order(ptSling[, i], na.last = NA)
    return(data.frame(cells= row.names(ptSling[idx,]),
                      continuous=ptSling[idx, i]))
  }))
}
