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

# Cardinality <- function(dataFrame) {
#     n <- dim(dataFrame)[2]
#     Df <- c()
#     for (i in 1:n) {
#         Df[i] <- length(unique(dataFrame[,i]))
#     }
#     return (Df)
# }

Statistic <- function(dataFrame, target, selected, subset, char=TRUE) {
  type <- 1
  if (char)
      type <- 0
  out <- 0
  ABCvec <- as.vector(as.matrix(dataFrame[1, c(target, selected, subset)]))
  ACvec <- as.vector(as.matrix(ABCvec[c(1, 3:length(ABCvec))]))
  BCvec <- as.vector(as.matrix(ABCvec[2:length(ABCvec)]))
  Cvec <- as.vector(as.matrix(ABCvec[3:length(ABCvec)]))
  tmp <- 1
  partialFrame <- dataFrame[, c(target, selected, subset)]
  while (tmp > 0) {
    ABClist <- Test(as.matrix(partialFrame),ABCvec,type)
    toBeDeletedItems <- ABClist[2:length(ABClist)]
    ABCvalue <- ABClist[1]
    ACvalue <- Test(as.matrix(dataFrame)[,c(target,subset)],ACvec,type)[1]
    BCvalue <- Test(as.matrix(dataFrame)[,c(selected,subset)],BCvec,type)[1]
    print("hi")
    print(class(as.matrix(dataFrame)[, subset]))
    Cvalue <- Test(as.matrix(dataFrame)[,subset],Cvec,type)[1]
    partialFrame <- partialFrame[toBeDeletedItems,]
    if(length(partialFrame)==length(c(target,selected,subset))) tmp <- 1
    else tmp <- dim(as.matrix(partialFrame))[1]
    if (tmp>0) {
      if(length(partialFrame)==length(c(target,selected,subset))) {
        ABCvec <- partialFrame
        tmp <- 0
      } else ABCvec <- partialFrame[1,]
          ACvec <- ABCvec[c(1,3:length(ABCvec))]
          BCvec <- ABCvec[2:length(ABCvec)]
          Cvec <- ABCvec[3:length(ABCvec)]
      }
    out <- out + 2*ABCvalue*log((ABCvalue*Cvalue)/(ACvalue*BCvalue))
  }
  TOL <- 0.001
  if(out < TOL) out <- 0
  return (out)
}

Pvalue <- function(tmp,a,b,c) {
    # tmp <- Example(1000)
    statistic <- Statistic(tmp,a,b,c)
    # df <- Cardinality(tmp)
    df <- Df(as.numeric(tmp$prop)[1:5],a,b,c)
    out <- pchisq(statistic, df, lower.tail = FALSE)
    return (out)
}