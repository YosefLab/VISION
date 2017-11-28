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
PseudotimeProjection = function(projection, root) {
  # find root cell tree vertex
  rootIdx = match(root, colnames(projection@pData))
  rootPos = projection@edgePos[rootIdx]
  if (rootPos >= 0.5) {
    rootVert = projection@edgeAssoc[2, rootIdx]
  }
  else {
    rootVert = projection@edgeAssoc[1, rootIdx]
  }
  paths = list()
  notVisited = c(1:length(projection@adjMat[1,]))



  findPaths = function(projection,
                       paths,
                       curVert,
                       prevVert,
                       curPath) {

    neighbors = c(which(projection@adjMat[curVert,]==1))
    toVisit = neighbors[which(neighbors != prevVert)]



    #check if end node
    if (length(toVisit) == 0) {
      #end of branch
      i = length(paths)+1
      paths[paste('path', i)] = list(curPath)
      return(paths)
    }



    while (length(toVisit) > 1) {
      # create new branch
      newPath = curPath

      nextVert = toVisit[1]
      #find cells on path
      if (nextVert > curVert) {
        edge = c(curVert, nextVert)
      }
      else {
        edge = c(nextVert, curVert)
      }

      cellIdx = which(apply(projection@edgeAssoc, 2, function(i) all(edge == i)))
      cells = colnames(projection@pData)[cellIdx]

      edgeVec = projection@vData[,edge[2]] - projection@vData[,edge[1]]

      cellPos = projection@edgePos[cellIdx]

      relProj = t(t(matrix(rep(edgeVec, length(cellPos)),
                           nrow=length(edgeVec),
                           ncol=length(cellPos)))*cellPos)


      #cellProj = relProj + projection@vData[,edge[1]]
      cellProj = projection@spatialPos[,cellIdx]

      #cellDist = curPath$distance + cellPos*sqrt(sum(edgeVec^2))
      if (is.null(dim(cellProj))) {
        cellDist = curPath$distance + sqrt(sum((cellProj-edgeVec)^2))
      } else {
        cellDist = curPath$distance + sqrt(colSums((cellProj-edgeVec)^2))
      }

      #add cells to path
      newPath$cells = c(curPath$cells, cells)
      newPath$cellProj = cbind(curPath$cellProj, cellProj)
      newPath$cellDist = c(curPath$cellDist, cellDist)
      newPath$distance = curPath$distance + sqrt(sum(edgeVec^2))


      paths = findPaths(projection,
                        paths,
                        nextVert,
                        curVert,
                        newPath)

      toVisit= toVisit[-1]


    }

    nextVert = toVisit[1]
    #find cells on path
    if (nextVert > curVert) {
      edge = c(curVert, nextVert)
    }
    else {
      edge = c(nextVert, curVert)
    }

    cellIdx = which(apply(projection@edgeAssoc, 2, function(i) all(edge == i)))
    cells = colnames(projection@pData)[cellIdx]

    edgeVec = projection@vData[,edge[2]] - projection@vData[,edge[1]]

    cellPos = projection@edgePos[cellIdx]

    relProj = t(t(matrix(rep(edgeVec, length(cellPos)),
                         nrow=length(edgeVec),
                         ncol=length(cellPos)))*cellPos)


    #cellProj = relProj + projection@vData[,edge[1]]
    cellProj = projection@spatialPos[,cellIdx]

    #cellDist = curPath$distance + cellPos*sqrt(sum(edgeVec^2))
    if (is.null(dim(cellProj))) {
      cellDist = curPath$distance + sqrt(sum((cellProj-edgeVec)^2))
    } else {
      cellDist = curPath$distance + sqrt(colSums((cellProj-edgeVec)^2))
    }


    #add cells to path
    curPath$cells = c(curPath$cells, cells)
    curPath$cellProj = cbind(curPath$cellProjs, cellProj)
    curPath$cellDist = c(curPath$cellDist, cellDist)
    curPath$distance = curPath$distance + sqrt(sum(edgeVec^2))


    paths = findPaths(projection,
                      paths,
                      nextVert,
                      curVert,
                      curPath)

    return(paths)
  }


  curPath = list(cells=c(), cellProj=c(), cellDist= c(), distance=0)
  paths = findPaths(projection,
                    paths,
                    rootVert,
                    rootVert,
                    curPath)

  maxDist = 0
  minDist = min(paths$`path 1`$cellDist)

  for (path in paths) {
    pathMax = max(path$cellDist)
    pathMin = min(path$cellDist)
    if (pathMax > maxDist) {
      maxDist = pathMax
    }
    if (pathMin < minDist) {
      minDist = pathMin
    }
  }
  distRange = maxDist - minDist

  pTimes = list()
  for (path in paths) {
    i = length(pTimes) + 1
    pts = (path$cellDist - minDist) / distRange
    df = data.frame(cell = path$cells, continuous=pts)
    rownames(df) = df$cell
    pTimes[paste('path',i)] = list(df)
  }
  return(pTimes)

}
#' Computes diffusion pseudotime.
#'
#' @param counts TreeProjection matrix
#' @param k # of neighbors for kernel width. Default sqrt(n).
#' @return list:
#' \itemize{
#'     \item path i: a list containing the cells and pseudotime coordinates
#'     for each path originating from root.
#' }
diffusionTime <- function(counts, root, k =NULL){
  nCells = NCOL(counts)
  rootIdx = match(root, colnames(counts))
  if (is.null(k)) {
    k = as.integer(sqrt(NCOL(counts))) + 1
  }
  distMat = sqdist(t(counts), t(counts))
  width = t(apply(distMat, 1, sort))[,k]

  sig = matrix(rep(width, nCells), nrow = nCells, ncol = nCells)
  sigCross = width %*% t(width)
  sigSum = sig + t(sig)
  kernel = (2*sigCross/sigSum)^0.5 * exp(-1/2*distMat/sigSum)
  Z = apply(kernel, 1, sum)
  Z = Z %*% t(Z)
  W = kernel / Z
  Zx = apply(W, 1, sum)
  Zy = apply(W, 2, sum)
  Z = Zx^-0.5 %*% t(Zy)^-0.5
  Tm = W*Z
  eig = eigen(Tm)
  rootPi= rep(0, nCells)
  rootPi[rootIdx] = 1

  dists = rowSums((eig$values[-1]/(1-eig$values[-1]))^2 * (c(eig$vectors[rootPi,-1])
                                                           - eig$vectors[,-1])^2)

}
