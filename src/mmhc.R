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
    cpc <- c()
    target <- 1
    iterate <- 2:dim(mat)[2]
    bads <- c()
    goods <- c()
    values <- c()
    maxP <- 0.0

    cpc2 <- c()

    while (target != 6) {

        for (i in iterate) {
            tmp <- mat[, c(target, i, cpc)]
            df <- Cardinality(tmp)
            pvalue <- Statistic(tmp, unique(tmp), df)
            if (pvalue > 0.05) {
                # bads <- c(bads, -i)
            } else {
                goods <- c(goods, i)
                # values <- c(values, pvalue)
            }
        }
        # u <- cbind(goods, values)
        # cpc <- c(cpc, which(max(u[,2]) == u)-dim(u)[1])
        cpc <- c(cpc, goods[1])
        iterate <- goods[goods != cpc]
        goods <- c()
        # bads <- c()
        maxP <- 0.0
        # target <- target
        if(length(iterate)==0) {
            cpc2 <- c(cpc2, cpc, 0)
            cpc <- c()
            target <- target + 1
            iterate <- 1:dim(mat)[2]
            iterate <- iterate[iterate != target]
        }
    }

    # goods <- c()
    # bads <- c()
    # iterate <- sample(1:dim(mat)[2])
    # target <- sample(iterate,1)
    # iterate <- 1:dim(mat)[2]
    # nots <- c(cpc, target)
    # iterate <- iterate[iterate != nots]
    # cpc <- cpc[cpc != target]
    # powSet <- 2 ^ cpc
    # powSet <- as.vector(powSet)

    # for (i in iterate) {
    #     tmp <- mat[, c(target, i, cpc)]
    #     df <- Cardinality(tmp)
    #     pvalue <- Statistic(tmp, unique(tmp), df)

    #     if (pvalue < 0.05)
    #         bads <- c(bads, i)
    #     else {
    #         goods <- c(goods, pvalue)
    #     }
    #     cpc <- c(cpc, pvalue)
    # }

    # mat <- as.matrix(Example(100))[, c(target, selected, subset)]
    # df <- Cardinality(mat)
    # out <- Statistic(mat, unique(mat), df)
    return (cpc2)
}