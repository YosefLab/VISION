library(Biobase)
library(SummarizedExperiment)

setwd("/Users/MattJones/FastProject_R/FastProject_R")

data <- as.matrix(read.table("LPS_and_Unstim_all.txt", sep="\t", header=TRUE, row.names=1))

expr <- ExpressionSet(data)
expr

se <- SummarizedExperiment(expr)
se