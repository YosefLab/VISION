% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-Signature.R
\name{sigConsistencyScores}
\alias{sigConsistencyScores}
\title{Evaluates the significance of each signature in each cluster}
\usage{
sigConsistencyScores(weights, sigScoresData, metaData, randomSigData, normExpr)
}
\arguments{
\item{weights}{output of computeKNNWeights}

\item{sigScoresData}{numeric matrix of signature scores
size is cells x signatures}

\item{metaData}{data.frame of meta-data for cells}

\item{randomSigData}{A list with two items:

  randomSigs: a list of signature score matrices.  Each list item
    represents permutation signatures generated for a specific size/balance,
    and is a numeric matrix of size cells X signatures

  sigAssignments: named factor vector assigning signatures to random background
    groups}

\item{normExpr}{NormData object used to generate random background sigs}
}
\value{
list:
\itemize{
    \item sigProbMatrix: the vector of consistency z-scores
    \item pVals: pvalues for the scores
    \item emp_pVals: pvalues for the scores
}
}
\description{
Evaluates the significance of each signature in each cluster
}
