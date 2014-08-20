library("Rcpp")
library("rbenchmark")
library("bnlearn")
library("igraph")
source('myExample.R')
sourceCpp('../src/mmpc.cpp')
sourceCpp('../src/score.cpp')

df <- student(1000)

MMHC <- function(df, alpha = 0.05) {
	columnNames <- colnames(df)
	mat <- data.matrix(df)
	PC <- MMPC(mat, alpha)
	adjMat <- BDeu(mat, PC, as.integer(dim(mat)[2]))
	colnames(adjMat) <- columnNames
	adjMat <- graph.adjacency(adjMat)

	return (adjMat)
}

print(MMHC(df))
rm(list = ls())

# T <- benchmark(mmhc(df, test = "x2", score = "bde"), MMHC(df), replications=1)