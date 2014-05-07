# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

if (!("Rcpp" %in% rownames(installed.packages())))
  install.packages("Rcpp")
if(!("sets" %in% rownames(installed.packages())))
  install.packages("sets")

require("Rcpp")
require("sets")
sourceCpp("mmhc.cpp")
source("mmhc_test.R")

Cardinality <- function(observedMatrix) {
    n <- dim(observedMatrix)[2]
    Df <- c()
    for (i in 1:n) {
        Df[i] <- length(unique(observedMatrix[,i]))
    }
    return (Df)
}

Run <- function(mat) {
    # mat <- as.matrix(Example(100))[, c(target, selected, subset)]
    df <- Cardinality(mat)
    out <- Statistic(mat, unique(mat), df)
    return (out)
}