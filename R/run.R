library("Rcpp")
library("bnlearn")
library("rbenchmark")
library("igraph")

source("myExample.R")

sourceCpp("../src//mmhc.cpp")

# df <- student(20000)

MMHC_R <- function(df) {
    columnNames <- colnames(df)
    K <- new(MMHC, df)
    K$mmpc()
    K$mmhc()
    adjMat <- K$adjMat()
    colnames(adjMat) <- columnNames
    adjMat <- graph.adjacency(adjMat)    
    return (adjMat)
}

MMPC_R <- function(df) {
    K <- new(MMHC, df)
    K$mmpc()
    return (K$pc())
}

# K <- new(MMHC, df)

# nobs <- c()
# timeM <- rep(0, 20)
# timeBN <- rep(0, 20)
# 
# for (j in 1:10) {
#     for (i in 1:20) {
#         nobs[i] <- 1000*i
#         df <- student(nobs[i])
#         K <- new(MMHC, df)
#         bm <- benchmark(mmpc(df), K$mmpc(),
#                         columns = c("test", "elapsed", "relative"),
#                         replications = 1)
#         timeM[i] <- timeM[i] + bm["elapsed"][[1]][1]
#         timeBN[i] <- timeBN[i] + bm["elapsed"][[1]][2]
#     }
#     print(j)
# }
# 
# timeM <- timeM/10
# timeBN <- timeBN/10
# 
# plot(nobs, timeM, type="l", col = "red", xlab = "number of observations", ylab = "time in s")
# lines(nobs, timeBN, col = "green")

nobs <- c()
timeM <- rep(0, 10)
timeBN <- rep(0, 10)

for (j in 1:10) {
    for (i in 1:10) {
        nobs[i] <- 1000*i
        df <- student(nobs[i])
        K <- new(MMHC, df)
        K$mmpc()
        bm <- benchmark(mmhc(df, score = "bde"), K$mmhc(),
                        columns = c("test", "elapsed", "relative"),
                        replications = 1)
        timeM[i] <- timeM[i] + bm["elapsed"][[1]][1]
        timeBN[i] <- timeBN[i] + bm["elapsed"][[1]][2]
    }
    print(j)
}

timeM <- timeM/10
timeBN <- timeBN/10

plot(nobs, timeBN, type="l", col = "red", xlab = "number of observations", ylab = "time in s")
lines(nobs, timeM, col = "green")

# bm <- benchmark(mmpc(df), K$mmpc(),
#                 columns = c("test", "elapsed", "relative"),
#                 replications = 1)