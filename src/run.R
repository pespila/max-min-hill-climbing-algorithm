library("Rcpp")
library("bnlearn")
library("rbenchmark")
library("igraph")
sourceCpp("mmhc.cpp")

df <- student(1000)

my <- function(df, alpha = 0.05, eta = as.integer(1)) {
	K <- new(MMHC, df, alpha, eta)
	K$mmpc()
	K$mmhc()
	adjMat <- graph.adjacency(K$graph)
	print(K$score)
	# plot(adjMat)
}

# T <- benchmark(mmhc(df), my(df), replications=1)