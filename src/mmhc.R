library("Rcpp")
library("bnlearn")
library("rbenchmark")
library("igraph")
sourceCpp('mmpc.cpp')
sourceCpp('score.cpp')

MMHC <- function(df, alpha = 0.05) {
	columnNames <- colnames(df)
	mat <- data.matrix(df)
	PC <- MMPC(mat, alpha)
	adjMat <- BDeu(mat, PC, as.integer(dim(mat)[2]))
	colnames(adjMat) <- columnNames
	adjMat <- graph.adjacency(adjMat)

	return (adjMat)
}

# T <- benchmark(mmhc(df, test = "x2", score = "bde"), MMHC(df), replications=1)