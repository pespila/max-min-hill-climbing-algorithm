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

card <<- c(2, 2, 2, 3, 2)
# Matrix <- as.matrix(Example(250, char = FALSE))

MaxMinHeuristic <- function(T, CPC, Matrix, maxNumberOfVariables, selectedBefore = 0, minimum = 1) {
	reject <- c(selectedBefore, minimum, 1)
	temporaryMinimum <- reject
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
				temporaryMinimum <- reject
				reject <- c(X, pvalue[1], pvalue[2])
			} else if (pvalue[1] == reject[2] && pvalue[2] > reject[3]) {
				reject <- c(X, reject[2], pvalue[2])
			}
		} else {
			accepted <- c(accepted, X)
		}
	}
	# out <- list("accepted" = accepted, "CPC" = reject, "tmpMin" = temporaryMinimum)
	out <- list()
	if (reject[1] == 0) {
		out <- list("accepted" = accepted)
	} else if (temporaryMinimum[1] != 0) {
		out <- list("accepted" = accepted, "CPC" = reject, "tmpMin" = temporaryMinimum)
	} else {
		out <- list("accepted" = accepted, "CPC" = reject)
	}
	return (out)
}


ForwardPhase <- function(T, Matrix) {
	temporaryMinimum <- 0
	maxNumberOfVariables <- 1:dim(Matrix)[2]
	maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == T)]
	CPCset <- MaxMinHeuristic(T, NULL, Matrix, maxNumberOfVariables)
	CPC <- as.integer(CPCset$CPC[1])
	crossOuts <- c(as.integer(CPCset$accepted), CPC)
	
	for (crossOut in crossOuts) {
		maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
	}

	repeat {
		if (length(maxNumberOfVariables) == 0) {
			break
		}

		if (length(CPC) == 1) {
			if (length(CPCset$tmpMin) == 0) {
				CPCset <- MaxMinHeuristic(T, CPC, Matrix, maxNumberOfVariables)
			} else {
				CPCset <- MaxMinHeuristic(T, CPC, Matrix, maxNumberOfVariables, CPCset$tmpMin[1], CPCset$tmpMin[2])
			}
			# CPCset <- MaxMinHeuristic(T, CPC, Matrix, maxNumberOfVariables)
			CPC <- c(CPC, as.integer(CPCset$CPC[1]))
			temporaryMinimum <- CPCset$CPC[2]
			crossOuts <- c(as.integer(CPCset$accepted), CPC)

			for (crossOut in crossOuts) {
				maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
			}
		} else {
			CPCiterationSet <- 2 ^ as.set(CPC)
			tmpCPCset <- 2 ^ as.set(CPC[1:(length(CPC)-1)])
			CPCiterationSet <- CPCiterationSet - tmpCPCset
			CPCiterationList <- as.list(CPCiterationSet)
			n <- 1
			rejectList <- list()
			reject <- c(0, 1)
			for (cpc in CPCiterationList) {
				if (length(CPCset$tmpMin) == 0 || CPCset$tmpMin[1] %in% CPC) {
					CPCset <- MaxMinHeuristic(T, as.numeric(cpc), Matrix, maxNumberOfVariables)
				} else {
					CPCset <- MaxMinHeuristic(T, as.numeric(cpc), Matrix, maxNumberOfVariables, CPCset$tmpMin[1], CPCset$tmpMin[2])
				}

				if (length(CPCset$CPC) != 0) {
					rejectList[[n]] <- CPCset$CPC
					n <- n + 1
				}
			}

			if (length(rejectList) != 0) {
				for (x in rejectList) {
					if (x[2] < reject[2]) {
						reject <- x
					}
				}
				CPC <- c(CPC, as.numeric(reject[1]))
				temporaryMinimum <- reject[2]
				crossOuts <- c(as.integer(CPCset$accepted), CPC)
				for (crossOut in crossOuts) {
					maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
				}
			} else {
				maxNumberOfVariables <- c()
			}
				# maxNumberOfVariables <- c()
		}

	}

	return (CPC)
}

# for (i in 1:5) {
# 	print(ForwardPhase(i, Matrix))
# }