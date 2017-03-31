#' Test file for generating plots of output from FastProject to check against other iterations of the program

projectionPlots <- function(projData) {
  
  for (p in projData@projections) {
    x <- p@pData[1,]
    y <- p@pData[2,]
    plot(x, y)
    title(paste0(projData@filter, ", ", p@name))
  }
  
}