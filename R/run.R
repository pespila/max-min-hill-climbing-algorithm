library("Rcpp")
library("bnlearn")
library("rbenchmark")
library("igraph")
sourceCpp("../src/mmhc.cpp")
source("myExample.R")

df <- student(1000)

my <- function(df) {
	K <- new(MMHC, df)
	K$mmpc()
	K$mmhc()
	adjMat <- graph.adjacency(K$adjMat())
	print(K$score())
	rm(K)
	# plot(adjMat)
}


# T <- benchmark(mmhc(df), my(df), replications=1)