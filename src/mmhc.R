# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

if (!("Rcpp" %in% rownames(installed.packages())))
  install.packages("Rcpp")
# if(!("sets" %in% rownames(installed.packages())))
#   install.packages("sets")

require("Rcpp")
require("bnlearn")
require("rbenchmark")
# require("sets")
sourceCpp("mmhc.cpp")
source("mmhc_test.R")

MaxMinHeuristic <- function(T, CPC, mat, card, iterate) {
    reject <- c(0, 0)
    accepted <- c()
    for (X in iterate) {
        print(c(X, T, CPC))
        statisticMatrix <- mat[, c(X, T, CPC)]
        pvalue <- Statistic(statisticMatrix, unique(statisticMatrix), card[c(X, T, CPC)])
        if (pvalue < 0.05) {
            if (pvalue > reject[2]) {
                reject <- c(X, pvalue)
            }
        } else {
            accepted <- c(accepted, X)
        }
    }
    reject <- c(reject, accepted)
    return (reject)
}

MMPCtmp <- function(mat, card, T) {
    CPC <- c()
    Fset <- c()
    iterate <- 1:dim(mat)[2]
    iterate <- iterate[!(iterate == T)]
    tmp <- c()

    repeat {
        tmp <- c()
        Fset <- MaxMinHeuristic(T, CPC[length(CPC)], mat, card, iterate)
        iterate <- iterate[!(iterate == Fset[1])]
        if (length(Fset) > 2) {
            tmp <- Fset[3:length(Fset)]
            iterate <- iterate[!(iterate == tmp)]
        }
        if (Fset[2] != 0) {
            CPC <- c(CPC, Fset[1])
        }
        if (length(tmp) == 0) {
            break
        }
    }
    return (CPC)
}

ForwardPhase <- function(mat, card, dimMat = dim(mat)[2], alpha = 0.05) {
    cpc <- c()
    target <- 1
    iterate <- 2:dim(mat)[2]
    bads <- c()
    goods <- c()
    values <- c()
    maxP <- 0.0
    n <- 1

    cpc2 <- matrix(0,dim(mat)[2],dim(mat)[2])

    while (target != dim(mat)[2]+1) {

        for (i in iterate) {
            tmp <- mat[, c(target, i, cpc[length(cpc)])]
            pvalue <- Statistic(tmp, unique(tmp), card[c(target, i, cpc[length(cpc)])])
            if (pvalue < 0.05) {
                goods <- c(goods, i)
                # print(goods)
            }
            # print(c(target, i, cpc[length(cpc)], i, pvalue))
        }
        cpc <- c(cpc, goods[1])
        iterate <- goods[goods != goods[1]]
        goods <- c()
        maxP <- 0.0
        if(length(iterate)==0) {
            for (i in (length(cpc)+1):dim(mat)[2]) {
                cpc <- c(cpc, 0)
            }
            cpc2[n,] <- cpc
            cpc <- c()
            n <- n + 1
            target <- target + 1
            iterate <- 1:dim(mat)[2]
            iterate <- iterate[iterate != target]
        }
    }
    return (cpc2)
}

BackwardPhase <- function(mat, cpc, card) {
    temp <- which(cpc[1,]!=0)
    tmp <- cpc[1,temp]



    for (i in 1:dim(cpc)[2]) {
        iterate <- 1:dim(mat)[2]
        stop <- 0
        for (j in 1:dim(cpc)[1]) {
            if (cpc[i, j] == 0) {
                next
            } else {
                iterate <- iterate[!(iterate == cpc[i, j])]
                for (T in iterate) {
                    for (k in 1:dim(cpc)[2]) {
                        for (l in 1:dim(cpc)[1]) {
                            if (cpc[k, l] == 0)  {
                                break
                            } else {
                                tmp <- mat[, c(cpc[i, j], T, cpc[k, l])]
                                pvalue <- Statistic(tmp, unique(tmp), card[c(cpc[i, j], T, cpc[k, l])])
                                if (pvalue < 0.05) {
                                    stop <- 1
                                    cpc[k, l] <- 0
                                } else {
                                }
                            }
                            if (stop) {
                                break
                            }
                        }
                        if (stop) {
                            break
                        }
                    }
                    if(stop) {
                        break
                    }
                }
            }
            if (stop) {
                break
            }
        }
    }
    return (cpc)
}

tempo <- as.matrix(Example(1000, char=FALSE))
card <- c(2,2,2,3,2)
for (i in 1:5) {
    print(MMPCtmp(tempo, card, i))
}
# FWP <- ForwardPhase(tempo, card)
# BWP <- BackwardPhase(tempo, FWP, card)
# test <- recursion(tempo, 1, 2:5)
# bm <- benchmark(Run(tempo), RunI(tempo), replications = 100, columns = c("elapsed", "relative", "test"))
# u<-c(2,2,2,3,2)
# a<-benchmark(Run(tempo),replications=1)
# v<-benchmark(Statistic(tempo[,c(1,2)], unique(tempo[,c(1,2)])),replications=29)
# x<-benchmark(Run(tempo),replications=1)
# x<-benchmark(Df(u),replications=1000)
# y<-benchmark(Cardinality(tempo),replications=1000)
# z<-benchmark(allC(1:10,1:10),replications=1000)