# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

# Checks whether 'Rcpp' is installed or not. If not it will be installed
if (!("Rcpp" %in% rownames(installed.packages())))
  install.packages("Rcpp")

# Checks whether 'sets' is installed or not. If not it will be installed
if(!("sets" %in% rownames(installed.packages())))
  install.packages("sets")

# Load the packages which are neede:
# - Rcpp: for compiling
# - bnlearn: will be removed
# - rbenchmark: will be removed
# - sets: to create and handle sets
# - mmhc.cpp: the C++-file which holds the 'fast' code blocks
# - mmhc_test.R: my Examples based on the book of Daphne Koller
require("Rcpp")
require("bnlearn")
require("rbenchmark")
require("sets")
sourceCpp("mmhc.cpp")
source("mmhc_test.R")

# Function MaxMinHeuristic which takes:
# - the target variable T for whose children and parents we are seeking for.
# - the current set of children and parents CPC
# - the maximum numbers of variables where we want to iterate over (maxNumberOfVariables)
# - the defaults: selectedBefore and minimum are 0 and 1 respectively. If there was a call of this
#   function before they could get set for statistical purposes.
# It returns a list of parents and children for a specific target variable, the variables where the
# nullhypothesis holds (they are going to be crossed out) and maybe the variables with the (second) highest association
# for which the nullhypothesis also could have been reject.

MaxMinHeuristic <- function(T, CPC, Matrix, maxNumberOfVariables, selectedBefore = 0, minimum = 1) { # MAXMINHEURISTIC
	reject <- c(selectedBefore, minimum, 1) # set which holds the values to compare (reject the nullhypothesis)
	temporaryMinimum <- reject
	accepted <- c() # the sets which holds the values where the nullhypothesis is accepted

	# we do not iterate over values in CPC because of computational costs
	for (crossOuts in CPC) { # FOR: crossing out
		
		maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOuts)]	
	
	} # FOR

	# Computing starts here
	for (X in maxNumberOfVariables) { # FOR: iteration over X
		
		setForSvalues <- c(X, T, CPC) # set the S values
		statisticMatrix <- Matrix[, setForSvalues] # take the specific columns of the matrix
		pvalue <- Statistic(statisticMatrix, unique(statisticMatrix), card[setForSvalues]) # compute the pvalue

		# statistical testing
		if (pvalue[1] < alpha) { # IF: reject nullhypothesis

			if (pvalue[1] < reject[2]) { # IF: test for the maximum over all rejected ones

				temporaryMinimum <- reject
				reject <- c(X, pvalue[1], pvalue[2])

			} else if (pvalue[1] == reject[2] && pvalue[2] > reject[3]) { # ELSE IF: if there is equality between two rejected ones take the G² value

				reject <- c(X, reject[2], pvalue[2])

			} # ELSE IF

		} else { #ELSE: accept nullhypothesis 

			accepted <- c(accepted, X)
		
		} # ELSE

	} # FOR

	# Set the output of this function
	out <- list()

	if (reject[1] == 0) { # IF: nothing to reject (in this case this does not happen)

		out <- list("accepted" = accepted)
	
	} else if (temporaryMinimum[1] != 0) { # ELSE IF: returns accepted (accept nullhypothesis), reject (CPC set) and the highest value after the maximum value
		
		out <- list("accepted" = accepted, "CPC" = reject, "tmpMin" = temporaryMinimum)
	
	} else { # ELSE: only return accepted and reject because there was no other value for which the nullhypothesis could be rejected

		out <- list("accepted" = accepted, "CPC" = reject)
	
	} # ELSE

	return (out)

} # MAXMINHEURISTIC

# Function ForwardPhase which takes the target variable T and the underlying matrix (later on a data frame)
# it returns a possible CPC set but it could have some false positive values

ForwardPhase <- function(T, Matrix) { # FORWARDPHASE

	temporaryMinimum <- 0 # will be set if there where more then one values for which the nullhypothesis could have been rejected
	maxNumberOfVariables <- 1:dim(Matrix)[2] # iteration array...
	maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == T)] # ...iterate over all values except the target (trivial case)
	CPCset <- MaxMinHeuristic(T, NULL, Matrix, maxNumberOfVariables) # the first CPC set where we start with the empty set
	CPC <- as.integer(CPCset$CPC[1]) # convert the CPC set as an array from the set (list) getting from one line above
	crossOuts <- c(as.integer(CPCset$accepted), CPC) # sets the variables where we do not iterate over again because they where accepted or rejected
	
	# set new iteration array
	for (crossOut in crossOuts) { # FOR

		maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
	
	} # FOR

	# main loop
	repeat { # REPEAT

		# break point: if there is nothing to iterate over the repeat breaks
		if (length(maxNumberOfVariables) == 0) { # IF

			break
		
		} # IF

		# If CPC has length one, we just want to take the last added value of CPC for calculation, else: take the subsets
		if (length(CPC) == 1) { # IF
			
			if (length(CPCset$tmpMin) == 0) {
				CPCset <- MaxMinHeuristic(T, CPC, Matrix, maxNumberOfVariables)
			} else {
				CPCset <- MaxMinHeuristic(T, CPC, Matrix, maxNumberOfVariables, CPCset$tmpMin[1], CPCset$tmpMin[2])
			}
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
		}

	} # REPEAT

	return (CPC)
}

# Function BackwardPhase which detects false positive elements of CPC and removes them
# Takes the current target T and the current CPC set of T. Returns the 'clean' CPC set

BackwardPhase <- function(T, CPC) { # BACKWARDPHASE

	if (length(CPC) == 1) { # IF

		return (CPC) # return the CPC set, if its length is 1. Trivial case...

	} else { # ELSE

		CPCset <- 2 ^ as.set(CPC) # power set of CPC
		CPCList <- as.list(CPCset) # make a list out of the power set for iteration and computing purpose
		remove <- c() # the array which should hold the values 'to be removed' from CPC

		# iterate over all values in CPC
		for (X in CPC) { # FOR

			# all 'combinations' from the CPCList (power set)
			for (cpc in CPCList) { # FOR

				setForSvalues <- c(X, T, as.numeric(cpc)) # S values for calculation of the pvalue
				statisticMatrix <- Matrix[, setForSvalues] # the underlying matrix for calculation
				pvalue <- Statistic(statisticMatrix, unique(statisticMatrix), card[setForSvalues]) # calculate the pvalue

				# if there exists a subset S of CPC for which a value X from CPC is removed, remove it and stop the for loop
				if (pvalue[1] > alpha && pvalue[1] != 1) { # IF

					CPC <- CPC[!(CPC == X)]
					break
				
				} # IF
			} # FOR
		} # FOR

		return (CPC)

	} # ELSE
} #BACKWARDPHASE

# Function MMPC. Takes the observed matrix (later on a data frame) and returns the 'real' set of parents and children for all
# possible target variables

MMPC <- function(Matrix) { # MMPC
	PC <- list()

	# iterate over all target variables
	for (T in 1:dim(Matrix)[2]) { # FOR

		CPC <- ForwardPhase(T, Matrix) # Runs the ForwardPhase
		CPC <- BackwardPhase(T, CPC) # Runs the BackwardPhase

		# calculates whether the observed target variable is a member of the temporary CPC set with target X out of CPC
		# if so then X from CPC is a parent or child of T, if not then it's not likely that X from CPC is a parent or child
		# of T and is removed.

		for (X in CPC) { # FOR

			tmp <- ForwardPhase(X, Matrix) # ForwardPhase for temporary CPC set
			tmp <- BackwardPhase(X, tmp) # BackwardPhase for temporary CPC set

			# test whether T is in temporary CPC set or not. If not -> remove corresponding X
			if (!(T %in% tmp)) { # IF
				
				CPC <- CPC[!(CPC == X)]
			
			} # IF
		}

		# Set PC
		PC[[T]] <- CPC

	} # FOR

	return (PC)
} # MMPC

# Actually not :)
Scoring <- function(PC) {
	scoreMatrix <- matrix(0, length(PC), length(PC))
	for (i in 1:length(PC)) {
		if (length(PC[[i]]) == 1) {
			scoreMatrix[i, PC[[i]][1]] <- 1
		}
	}
	return (scoreMatrix)
}

# library("igraph")
# mat <- matrix(0,5,5)
# mat[1,4]<-1
# mat[2,3]<-1
# mat[2,4]<-1
# mat[4,5]<-1
# aj <- graph.adjacency(mat)
# plot(aj)
# mat

# bench <- benchmark(MMPC(Matrix), mmpc(Example(250,char=FALSE)), replications=1)