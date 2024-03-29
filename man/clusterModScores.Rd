% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-Module.R
\name{clusterModScores}
\alias{clusterModScores}
\title{Compute Ranksums Test, for all factor meta data.  One level vs all others}
\usage{
clusterModScores(object, variables = "All")
}
\arguments{
\item{object}{the VISION object}

\item{variables}{which columns of the meta-data to use for comparisons}
}
\value{
the VISION object with the @ClusterComparisons modules slot populated
}
\description{
Compute Ranksums Test, for all factor meta data.  One level vs all others
}
