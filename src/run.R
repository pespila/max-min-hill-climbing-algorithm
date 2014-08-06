library("Rcpp")
library("igraph")
library("rbenchmark")
sourceCpp("code.cpp")
source("../R/mmhc.R")

df <- Example(5000, char = FALSE)

mmhc <- function(df) {
	dm <- data.matrix(df)
	N <- as.integer(dim(dm)[2])
	M <- MMPC(dm, 0.05, N)
	AdjMat <- graph.adjacency(M)
	return (AdjMat)
}

T <- benchmark(mmhc(df), C_MMHC(df), replications=1)
# K <- new(Properties, dm)
# K <- MMPC(df, 0.05)
# M <- MMPC(K, 0.05)

# MMHC <- function(df) {
# 	dm <- data.matrix(dm)

# }