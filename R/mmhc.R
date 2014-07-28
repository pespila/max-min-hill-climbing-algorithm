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
source("mmhc_test.R")
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
MaxMinHeuristic <- function(T, CPC, mat, maxNumberOfVariables, alpha, selectedBefore = 0, minimum = 1) { # MAXMINHEURISTIC
	reject <- c(selectedBefore, minimum, 1) # set which holds the values to compare (reject the nullhypothesis)
	temporaryMinimum <- reject
	accepted <- c() # the sets which holds the values where the nullhypothesis is not rejected

	# we do not iterate over values in CPC because of computational costs
	for (crossOuts in CPC) { # FOR: crossing out
		
		maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOuts)]	
	
	} # FOR

	# Computing starts here
	for (X in maxNumberOfVariables) { # FOR: iteration over X
		
		setForSvalues <- c(X, T, CPC) # set the S values
		statisticMatrix <- mat[, setForSvalues] # take the specific columns of the matrix
		pvalue <- MySvalue(statisticMatrix) # compute the pvalue

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
		
		} # IF/ELSE

	} # FOR

	# Set the output of this function
	out <- list()

	if (reject[1] == 0) { # IF: nothing to reject (normally this does not happen)

		out <- list("accepted" = accepted)
	
	} else if (temporaryMinimum[1] != 0) { # ELSE IF: returns accepted (accept nullhypothesis), reject (CPC set) and the highest value after the maximum value
		
		out <- list("accepted" = accepted, "CPC" = reject, "tmpMin" = temporaryMinimum)
	
	} else { # ELSE: only return accepted and reject because there was no other value for which the nullhypothesis could be rejected

		out <- list("accepted" = accepted, "CPC" = reject)
	
	} # ELSE

	return (out)

} # MAXMINHEURISTIC

# Function ForwardPhase which takes the target variable T and the underlying matrix (later on a data frame)
# it returns a possible CPC set which could have some false positive values
ForwardPhase <- function(T, mat, alpha) { # FORWARDPHASE

	tmp <- list()
	CPC <- UpdateCPC(tmp, 0)
	maxNumberOfVariables <- 1:dim(mat)[2] # iteration array...
	maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == T)] # ...iterate over all values except the target (trivial case)
	CPCset <- MaxMinHeuristic(T, NULL, mat, maxNumberOfVariables, alpha) # the first CPC set where we start with the empty set
	if (length(CPCset) == 1) {
		maxNumberOfVariables <- NULL
	} else {
		CPC <- UpdateCPC(CPC, as.integer(CPCset$CPC[1]))
		crossOuts <- c(as.integer(CPCset$accepted), CPC[[length(CPC)]]) # sets the variables where we do not iterate over again because they where accepted or rejected #!!!!!!BEFOR CPC
		# set new iteration array
		for (crossOut in crossOuts) { # FOR

			maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
		
		} # FOR
	}

	# main loop
	repeat { # REPEAT

		# break point: if there is nothing to iterate over the repeat breaks
		if (length(maxNumberOfVariables) == 0) { # IF

			break
		
		} # IF

		# If CPC has length three, we just want to take the last added value of CPC for calculation, else: take the subsets
		if (length(CPC) == 3) { # IF !!!!!!!!!!!!!Before == 1

			# From chapter 6: decides if there was a minimum before which could also be taken
			if (length(CPCset) != 3) { # IF

				CPCset <- MaxMinHeuristic(T, CPC[[3]], mat, maxNumberOfVariables, alpha)
			
			} else { # ELSE

				CPCset <- MaxMinHeuristic(T, CPC[[3]], mat, maxNumberOfVariables, alpha, CPCset$tmpMin[1], CPCset$tmpMin[2])
			
			} # IF/ELSE
			
			# Set the CPC set, a possible temporary minimum and the values which gonna be crossed out for iteration

			if (length(CPCset$CPC) != 0)
				CPC <- UpdateCPC(CPC, as.integer(CPCset$CPC[1]))

			crossOuts <- c(as.integer(CPCset$accepted), CPC[[length(CPC)]])

			# Cross specific values out
			for (crossOut in crossOuts) { # FOR

				maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
			
			} # FOR

		} else { # ELSE

			# According to chapter 6: just iterate over the subsets where the new node (random variable) is element of.
			# First take the power set of CPC, then a second power set without the new added value (in last iteration), subtract
			# the temporary power set from the original one and make a list out of it (for computations)
			CPCiterationList <- CPC[CPC[[1]]:length(CPC)]

			n <- 1 # helps for iteration

			rejectList <- list() # list with the rejected nullhypothesises
			reject <- c(0, 1) # single rejecting array, needed for the temporary minimum of this iteration (if there exists one)
			
			# Go over all subset S (namely cpc) in CPC (CPCiterationList) and call the MaxMinHeuristic depending if there was another
			# minimum before or not.
			for (cpc in CPCiterationList) { # FOR

				if ( length(CPCset) != 3 || ( length(CPCset) == 3 && CPCset$tmpMin[1] %in% CPC[[length(CPC)]] ) ) { # IF

					CPCset <- MaxMinHeuristic(T, as.numeric(cpc), mat, maxNumberOfVariables, alpha)
				
				} else { # IF

					CPCset <- MaxMinHeuristic(T, as.numeric(cpc), mat, maxNumberOfVariables, alpha, CPCset$tmpMin[1], CPCset$tmpMin[2])
				
				} # IF/ELSE

				# If CPC is not empty then there was a variable to reject the nullhypothesis with
				if (length(CPCset$CPC) != 0) { # IF

					rejectList[[n]] <- CPCset$CPC
					n <- n + 1

				} # IF


			} # FOR

			# Proof whether there was a new variable for which the nullhypothesis could be rejected. Also test if there were two or more variables
			# for which the nullhypothesis could be rejected. The 'maximum' (in this computation it's actually a minimum) without the real maximum/minimum
			# is taken for the temporary minimum set.
			if (length(rejectList) != 0) { # IF

				# Testing for maximum pvalue (actually a minimum is computed because I took the positive and not the negative pvalue)
				for (x in rejectList) { # FOR

					# minimum testing
					if (x[2] < reject[2]) { # IF

						reject <- x
					
					} # IF
				
				} # FOR

				# Set CPC, temporary minimum and cross out values
				CPC <- UpdateCPC(CPC, as.integer(reject[1]))

				crossOuts <- c(as.integer(CPCset$accepted), CPC[[length(CPC)]])

				# cross specific values out
				for (crossOut in crossOuts) { # FOR

					maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == crossOut)]
				
				}
			
			} else { # ELSE

				maxNumberOfVariables <- c()

			} # IF/ELSE

		} # IF/ELSE

	} # REPEAT

	return (CPC)
}

# Function BackwardPhase which detects false positive elements of CPC and removes them
# Takes the current target T and the current CPC set of T. Returns the 'clean' CPC set
BackwardPhase <- function(T, CPCset, mat, alpha) { # BACKWARDPHASE

	CPC <- CPCset[[length(CPCset)]]
	CPCList <- CPCset[2:length(CPCset)]

	if (length(CPC) == 1) { # IF

		return (CPC) # return the CPC set, if its length is 1. Trivial case...

	} else { # ELSE

		remove <- c() # the array which should hold the values 'to be removed' from CPC

		# iterate over all values in CPC
		for (X in CPC) { # FOR

			# all 'combinations' from the CPCList (power set)
			for (cpc in CPCList) { # FOR

				setForSvalues <- c(X, T, as.numeric(cpc)) # S values for calculation of the pvalue
				statisticMatrix <- mat[, setForSvalues]
				pvalue <- MySvalue(statisticMatrix)

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
R_MMPC <- function(mat, alpha) { # MMPC
	PC <- list()

	# iterate over all target variables
	for (T in 1:dim(mat)[2]) { # FOR

		CPC <- ForwardPhase(T, mat, alpha) # Runs the ForwardPhase
		if (length(CPC) == 2) {
			PC[[T]] <- NULL
			next
		}
		CPC <- BackwardPhase(T, CPC, mat, alpha) # Runs the BackwardPhase
		# calculates whether the observed target variable is a member of the temporary CPC set with target X out of CPC
		# if so then X from CPC is a parent or child of T, if not then it's not likely that X from CPC is a parent or child
		# of T and is removed.

		for (X in CPC) { # FOR

			tmp <- ForwardPhase(X, mat, alpha) # ForwardPhase for temporary CPC set
			tmp <- BackwardPhase(X, tmp, mat, alpha) # BackwardPhase for temporary CPC set

			# test whether T is in temporary CPC set or not. If not -> remove corresponding X
			if (!(T %in% tmp)) { # IF
				
				CPC <- CPC[!(CPC == X)]
			
			} # IF
		}

		# Set PC
		PC[[T]] <- CPC

	} # FOR

	AdjMat <- matrix(0, dim(mat)[2], dim(mat)[2])
	for (i in 1:length(PC)) {
		for (element in PC[[i]]) {
			AdjMat[i, element] <- 1
		}
	}
	# AdjMat <- graph.adjacency(AdjMat)
	# return (PC)
	return (AdjMat)
} # MMPC

getQ <- function(mat, pc = NULL) {
	out <- 1
	
	for (i in pc) {
		out <- out * getR(mat[, i])
	}

	return (out)
}

getEta <- function(mat, std = TRUE) {
	if (std) {
		return (1)
	} else {
		avNoOfValPerVar <- 0

		for (i in 1:dim(mat)[2]) {
			avNoOfValPerVar <- avNoOfValPerVar + getR(mat[, i])
		
		}
		
		avNoOfValPerVar <- avNoOfValPerVar/dim(mat)[2]
		
		return (avNoOfValPerVar)
	}
}

getNijk <- function(mat, pc, wij, Xi, k) {
	if (length(pc) == 0) {
		return (length(which(mat[, Xi] == k)))
	} else {
		count <- 0
		K <- which(mat[, Xi] == k)
		for (i in K) {
			if (identical(as.vector(mat[i, pc]), as.vector(wij))) {
				count <- count + 1
			}
		}
		return (count)
	}
}

getNij <- function(mat, pc, wij, Xi, r) {
	sum <- 0
	for (k in 1:r) {
		sum <- sum + getNijk(mat, pc, wij, Xi, k)
	}
	return (sum)
}

getW <- function(mat, pc) {
	return (unique(mat[, pc]))
}

MatrixScore <- function(mat, PC) {
	finalScore <- c()
	n <- dim(mat)[2]
	nijk <- 0
	nij <- 0


	for (i in 1:5) {
		jSum <- c()
		kSum <- c()
		pc <- which(PC[, i] != 0)
		q <- getQ(mat, pc)
		r <- getR(mat[, i])
		eta <- getEta(mat, FALSE)
		wi <- getW(mat, pc)
		x <- eta / q
		y <- eta / ( r * q )
		if (length(pc) <= 1) {
			q <- length(wi)
		} else {
			q <- dim(wi)[1]
		}

		for (j in 1:q) {

			if (length(pc) <= 1) {
				nij <- getNij(mat, pc, wi[j], i, r)
			} else {
				nij <- getNij(mat, pc, wi[j, ], i, r)
			}

			for (k in 1:r) {

				if (length(pc) <= 1) {
					nijk <- getNijk(mat, pc, wi[j], i, k)
				} else {
					nijk <- getNijk(mat, pc, wi[j, ], i, k)
				}

				kSum[k] <- lgamma ( nijk + y ) - lgamma ( y )

			}

			jSum[j] <- lgamma ( x ) - lgamma ( nij + x ) + sum(kSum)

		}
		finalScore[i] <- sum(jSum)
	}

	return (sum(finalScore))
}

SetEdge <- function(mat, PC) {
	n <- length(PC)
	adjMat <- matrix(0, n, n)
	before <- MatrixScore(mat, adjMat)
	after <- 0
	for (i in 1:length(PC)) {
		for (j in 1:length(PC[[i]])) {
			if (adjMat[i, PC[[i]][j]] == 0 && adjMat[PC[[i]][j], i] == 0) {
				adjMat[i, PC[[i]][j]] = 1
				after <- MatrixScore(mat, adjMat)

				if (after > before) {
					before <- after
				} else {
					adjMat[i, PC[[i]][j]] = 0
				}
			}
		}
	}
	return (adjMat)
}

ReversingOrDeleting <- function(mat, PC, adjMat) {
	before <- MatrixScore(mat, adjMat)
	for (i in 1:length(PC)) {
		for (j in 1:length(PC[i])) {

			rnd <- sample(c(0, 1), 1, prob = c(0.3, 0.7))

			if (adjMat[i, PC[[i]][j]] == 0 && adjMat[PC[[i]][j], i] == 0) {
				adjMat[i, PC[[i]][j]] = 1
				after <- MatrixScore(mat, adjMat)

				if (after > before) {
					before <- after
				} else {
					adjMat[i, PC[[i]][j]] = 0
				}
			} else if (adjMat[PC[[i]][j], i] == 1 && rnd) {
				adjMat[i, PC[[i]][j]] = 1
				adjMat[PC[[i]][j], i] = 0
				after <- MatrixScore(mat, adjMat)

				if (after > before) {
					before <- after
				} else {
					adjMat[i, PC[[i]][j]] = 0
					adjMat[PC[[i]][j], i] = 1
				}
			} else if (adjMat[i, PC[[i]][j]] == 1 && !rnd) {
				adjMat[i, PC[[i]][j]] = 0
				after <- MatrixScore(mat, adjMat)

				if (after > before) {
					before <- after
				} else {
					adjMat[i, PC[[i]][j]] = 1
				}
			}

		}
	}

	return (adjMat)
}

BDeuR <- function(mat, PC) {
	adjMat <- SetEdge(mat, PC)

	for (i in 1:5) {
		adjMat <- ReversingOrDeleting(mat, PC, adjMat)
	}

	print(MatrixScore(mat, adjMat))

	return (adjMat)
}

R_MMHC <- function(df) {
	columnNames <- colnames(df)
	mat <- data.matrix(df)
	PC <- R_MMPC (mat, 0.05)
	print(PC)
	adjMat <- BDeuR(mat, PC)
	colnames(adjMat) <- columnNames
	adjMat <- graph.adjacency(adjMat)

	return (adjMat)
}

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