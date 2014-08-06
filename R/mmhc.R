# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

# Checks whether 'Rcpp' is installed or not. If not it will be installed
if (!("Rcpp" %in% rownames(installed.packages())))
  install.packages("Rcpp")

# Load the packages which are neede:
# - Rcpp: for compiling
# - bnlearn: will be removed
# - rbenchmark: will be removed
# - mmhc.cpp: the C++-file which holds the 'fast' code blocks
# - mmhc_test.R: my Example's based on the book of Daphne Koller
require("Rcpp")
require("bnlearn")
require("rbenchmark")
require("igraph")
require("RcppArmadillo")
# sourceCpp("./../src/mmhc.cpp")
sourceCpp("./../src/mmpc.cpp")
sourceCpp("./../src/score.cpp")
# source("./mmhc_test.R")
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

options(warn=-1)

# Function MaxMinHeuristic which takes:
# - the target variable T for whose children and parents we are seeking for.
# - the current set of children and parents CPC
# - the maximum numbers of variables where we want to iterate over (maxNumberOfVariables)
# - the defaults: selectedBefore and minimum are 0 and 1 respectively. If there was a call of this
#   function before they could get set for statistical purposes.
# It returns a list of parents and children for a specific target variable, the variables where the
# nullhypothesis holds (they are going to be crossed out) and maybe the variables with the (second) highest association
# for which the nullhypothesis also could have been reject.


# Function ForwardPhase which takes the target variable T and the underlying matrix (later on a data frame)
# it returns a possible CPC set which could have some false positive values


# Function BackwardPhase which detects false positive elements of CPC and removes them
# Takes the current target T and the current CPC set of T. Returns the 'clean' CPC set


# Function MMPC. Takes the observed matrix (later on a data frame) and returns the 'real' set of parents and children for all
# possible target variables


C_MMHC <- function(df) {
	columnNames <- colnames(df)
	mat <- data.matrix(df)
	card <- Cardinality(mat)
	PC <- C_MMPC(mat, card, 0.05)
	adjMat <- BDeu(mat, card, PC, as.integer(dim(mat)[2]))
	colnames(adjMat) <- columnNames
	adjMat <- graph.adjacency(adjMat)

	return (adjMat)
}