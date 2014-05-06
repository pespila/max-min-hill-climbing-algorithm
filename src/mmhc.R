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

# MMPC.DegreesOfFreedom <- function(Df,target,selected,subset) {
#     out <- (Df[target]-1)*(Df[selected]-1)
#     if(length(subset)==1) {
#         out <- out*Df[subset]
#     } else {
#         for (i in 1:length(subset)) {
#             out <- out*Df[subset[i]]
#         }
#     }
#     return (out)
# }

MMPC.Statistic <- function(observedMatrix,target,selected,subset) {
    out <- 0
    ABCvec <- observedMatrix[1,c(target,selected,subset)]
    ACvec <- ABCvec[c(1,3:length(ABCvec))]
    BCvec <- ABCvec[2:length(ABCvec)]
    Cvec <- ABCvec[3:length(ABCvec)]
    tmp <- 1
    partialMatrix <- observedMatrix[,c(target,selected,subset)]
    type <- 1
    if(class(observedMatrix[1,1])=="character") type<-0
    while(tmp>0) {
        ABClist <- Test(partialMatrix,ABCvec,type)
        toBeDeletedItems <- ABClist[2:length(ABClist)]
        ABCvalue <- ABClist[1]
        ACvalue <- Test(observedMatrix[,c(target,subset)],ACvec,type)[1]
        BCvalue <- Test(observedMatrix[,c(selected,subset)],BCvec,type)[1]
        Cvalue <- Test(as.matrix(observedMatrix[,subset]),Cvec,type)[1]
        partialMatrix <- partialMatrix[toBeDeletedItems,]
        if(length(partialMatrix)==length(c(target,selected,subset))) tmp <- 1
        else tmp <- dim(as.matrix(partialMatrix))[1]
        if(tmp>0) {
            if(length(partialMatrix)==length(c(target,selected,subset))) {
                    ABCvec <- partialMatrix
                    tmp <- 0
                }
            else ABCvec <- partialMatrix[1,]
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
    statistic <- Statistic(as.matrix(tmp[,1:(dim(tmp)[2]-1)]),a,b,c)
    # df <- Cardinality(tmp)
    df <- Df(as.numeric(tmp$prop)[c(a,b,c)])
    out <- pchisq(statistic, df, lower.tail = FALSE)
    return (out)
}