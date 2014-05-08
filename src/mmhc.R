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

Run <- function(mat) {
    cpc <- c()
    target <- 1
    iterate <- 2:dim(mat)[2]
    bads <- c()
    goods <- c()
    values <- c()
    maxP <- 0.0

    cpc2 <- c()

    while (target != dim(mat)[2]) {

        for (i in iterate) {
            tmp <- mat[, c(target, i, cpc)]
            pvalue <- Statistic(tmp, unique(tmp))
            if (pvalue < 0.05)
                goods <- c(goods, i)
        }
        cpc <- c(cpc, goods[1])
        iterate <- goods[goods != goods[1]]
        goods <- c()
        maxP <- 0.0
        if(length(iterate)==0) {
            cpc2 <- c(cpc2, cpc, 0)
            cpc <- c()
            target <- target + 1
            iterate <- 1:dim(mat)[2]
            iterate <- iterate[iterate != target]
        }
    }
    return (cpc2)
}

tempo <- as.matrix(Example(100, char=FALSE))

bench <- function(tmp) {
    Statistic(tmp, unique(tmp))
}

a<-benchmark(Run(tempo),replications=1)
x<-benchmark(bench(tempo),replications=29)
y<-benchmark(Cardinality(tempo),replications=1000)