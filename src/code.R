# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

if (!("Rcpp" %in% rownames(installed.packages())))
  install.packages("Rcpp")
if(!("sets" %in% rownames(installed.packages())))
  install.packages("sets")

require("Rcpp")
require("bnlearn")
require("rbenchmark")
require("sets")
sourceCpp("mmhc.cpp")
source("mmhc_test.R")

cardinality <<- c(2, 2, 2, 3, 2)
Matrix <- as.matrix(Example(250, char = FALSE))

MaxMinHeuristic <- function(T, CPC, Matrix, maxNumberOfVariables) {
	reject <- c(0, 1, 1)
	accepted <- c()
	for (crossOuts in CPC) {
		maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOuts)]	
	}
	for (X in maxNumberOfVariables) {
		setForSvalues <- c(X, T, CPC)
		statisticMatrix <- Matrix[, setForSvalues]
		pvalue <- Statistic(statisticMatrix, unique(statisticMatrix), card[setForSvalues])
		if (pvalue[1] < 0.05) {
			if (pvalue[1] < reject[2]) {
				reject <- c(X, pvalue[1], pvalue[2])
			}
			if (pvalue[1] == reject[2] && pvalue[2] > reject[3]) {
				reject <- c(X, reject[2], pvalue[2])
			}
		} else {
			accepted <- c(accepted, X)
		}
	}
	out <- list()
	if (reject[1] == 0) {
		out <- list("accepted" = accepted)
	} else {
		out <- list("accepted" = accepted, "CPC" = reject)
	}
	return (out)
}


ForwardPhase <- function(T, Matrix) {
	maxNumberOfVariables <- 1:dim(Matrix)[2]
	maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == T)]
	CPCset <- MaxMinHeuristic(T, NULL, Matrix, maxNumberOfVariables)
	CPC <- as.integer(CPCset$CPC[1])
	crossOuts <- c(as.integer(CPCset$accepted), CPC)
	
	for (crossOut in crossOuts) {
		maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
	}

	return (maxNumberOfVariables)
}

# print(MaxMinHeuristic(1, NULL, Matrix))
for (i in 1:5) {
	print(ForwardPhase(i, Matrix))
}