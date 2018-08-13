#' VISION
#'
#' VISION analyzes expression data and produces a dynamic report in which various two dimensional
#' projections of the data can be explored. Annotated gene sets ("signatures") are incorporated so projection features
#' can be understood in relation to the biological processes present in the data. In addition to standard dimensionality
#' reduction methods, VISION incorporates fitting a high dimensional principle tree to allow exploring continuous
#' processes. VISION supports the analysis of massive datasets - theoretically up to a million single cells at any
#' sequencing depth, and an interactive visualization platforms that enables subset analysis based on two dimensional
#' representations.
#'
#' @docType package
#' @import Rcpp
#' @importFrom Rcpp evalCpp sourceCpp
#' @import methods
#' @import jug
#' @import loe
#' @import Matrix
#' @useDynLib VISION
#' @name VISION
NULL
