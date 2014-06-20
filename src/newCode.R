# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

# Checks whether 'Rcpp' is installed or not. If not it will be installed
if (!("Rcpp" %in% rownames(installed.packages())))
  install.packages("Rcpp")

# Checks whether 'sets' is installed or not. If not it will be installed
# if(!("sets" %in% rownames(installed.packages())))
#   install.packages("sets")

# Load the packages which are neede:
# - Rcpp: for compiling
# - bnlearn: will be removed
# - rbenchmark: will be removed
# - sets: to create and handle sets
# - mmhc.cpp: the C++-file which holds the 'fast' code blocks
# - mmhc_test.R: my Example's based on the book of Daphne Koller
require("Rcpp")
require("RcppArmadillo")
require("bnlearn")
require("rbenchmark")
require("igraph")
# require("sets")
sourceCpp("newMMPC.cpp")
source("mmhc_test.R")
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
tmp <<- Example(1000, char=FALSE)
options(warn=-1)
# MatrixT <<- t(Matrixy)

# Function MaxMinHeuristic which takes:
# - the target variable T for whose children and parents we are seeking for.
# - the current set of children and parents CPC
# - the maximum numbers of variables where we want to iterate over (maxNumberOfVariables)
# - the defaults: selectedBefore and minimum are 0 and 1 respectively. If there was a call of this
#   function before they could get set for statistical purposes.
# It returns a list of parents and children for a specific target variable, the variables where the
# nullhypothesis holds (they are going to be crossed out) and maybe the variables with the (second) highest association
# for which the nullhypothesis also could have been reject.
MaxMinHeuristic <- function(T, CPC, mat, maxNumberOfVariables, selectedBefore = 0, minimum = 1) { # MAXMINHEURISTIC
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
ForwardPhase <- function(T, mat) { # FORWARDPHASE

	tmp <- list()
	CPC <- UpdateCPC(tmp, 0)
	maxNumberOfVariables <- 1:dim(mat)[2] # iteration array...
	maxNumberOfVariables <- maxNumberOfVariables[!(maxNumberOfVariables == T)] # ...iterate over all values except the target (trivial case)
	CPCset <- MaxMinHeuristic(T, NULL, mat, maxNumberOfVariables) # the first CPC set where we start with the empty set
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

				CPCset <- MaxMinHeuristic(T, CPC[[3]], mat, maxNumberOfVariables)
			
			} else { # ELSE

				CPCset <- MaxMinHeuristic(T, CPC[[3]], mat, maxNumberOfVariables, CPCset$tmpMin[1], CPCset$tmpMin[2])
			
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

					CPCset <- MaxMinHeuristic(T, as.numeric(cpc), mat, maxNumberOfVariables)
				
				} else { # IF

					CPCset <- MaxMinHeuristic(T, as.numeric(cpc), mat, maxNumberOfVariables, CPCset$tmpMin[1], CPCset$tmpMin[2])
				
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
BackwardPhase <- function(T, CPCset, mat) { # BACKWARDPHASE

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
MMPC <- function(mat) { # MMPC
	PC <- list()

	# iterate over all target variables
	for (T in 1:dim(mat)[2]) { # FOR

		CPC <- ForwardPhase(T, mat) # Runs the ForwardPhase
		if (length(CPC) == 2) {
			PC[[T]] <- NULL
			next
		}
		CPC <- BackwardPhase(T, CPC, mat) # Runs the BackwardPhase
		# calculates whether the observed target variable is a member of the temporary CPC set with target X out of CPC
		# if so then X from CPC is a parent or child of T, if not then it's not likely that X from CPC is a parent or child
		# of T and is removed.

		for (X in CPC) { # FOR

			tmp <- ForwardPhase(X, mat) # ForwardPhase for temporary CPC set
			tmp <- BackwardPhase(X, tmp, mat) # BackwardPhase for temporary CPC set

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

	return (AdjMat)
} # MMPC

getR <- function(vec) {
	return (length(unique(vec)))
}

getDim <- function(mat) {
	return (dim(mat)[2])
}

getQ <- function(mat, pc = NULL) {
	if (length(pc) == 0) {
		return (1)
	} else if (length(pc) == 1) {
		return (length(unique(mat[, pc])))
	} else {
		return (dim(unique(mat[, pc]))[1])
	}
}

getEta <- function(mat, std = TRUE) {
	if (std) {
		return (1)
	} else {
		hDim <- 1:getDim(mat)
		avNoOfValPerVar <- 0
		for (i in hDim) {
			avNoOfValPerVar <- avNoOfValPerVar + getR(mat[, i])
		}
		avNoOfValPerVar <- avNoOfValPerVar/getDim(mat)
		return (avNoOfValPerVar)
	}
}

getNijk <- function(mat, pc, wij, Xi, j, k) {
	if (length(pc) == 0) {
		return (length(which(mat[, Xi] == k)))
	} else {
		count <- 0
		K <- which(mat[, Xi] == k)
		L <- which(mat[, pc] == wij)
		count <- length(which(K == L))
		# for (i in K) {
		# 	tmp <- mat[i, pc]
		# 	if (identical(tmp, wij))
		# 		count <- count + 1
		# }
		return (count)
	}
}

getNij <- function(mat, pc, wij, Xi, j, r) {
	sum <- 0
	for (k in 1:r) {
		sum <- sum + getNijk(mat, pc, wij, Xi, j, k)
	}
	return (sum)
}

getW <- function(mat, pc) {
	return (unique(mat[, pc]))
}

nuScore <- function(mat, PC) {
	scores <- c()
	n <- getDim(mat)
	nijk <- 0
	nij <- 0
	iSum <- 0

	for (i in 1:n) {
		
		pc <- PC[[i]]
		q <- getQ(mat, pc)
		r <- getR(mat[, i])
		eta <- getEta(mat)
		wi <- getW(mat, pc)
		x <- eta / q
		y <- eta / ( r * q )
		jSum <- 0

		for (j in 1:q) {

			nij <- getNij(mat, pc, wi[j], i, j, r)
			kSum <- 0

			for (k in 1:r) {

				nijk <- getNijk(mat, pc, wi[j], i, j, k)
				kSum <- kSum + lgamma ( nijk + y ) - lgamma ( y )

			}

			jSum <- jSum + lgamma ( x ) - lgamma ( nij + x ) + kSum

			# print(jSum)
		}

		scores[i] <- jSum
		# iSum <- iSum + jSum

	}

	return (sum(scores))
}

Scoring <- function(mat, PC) {
	n <- dim(PC)[1]
	scoreList <- initEmptyList(n)
	tmpScoreList <- scoreList
	noChange <- 0
	scores <- nuScore(mat, scoreList)


	repeat {
		if (noChange == 40)
			break

		sampling <- sample(1:5, 2)

		if (!(sampling[2] %in% tmpScoreList[[sampling[1]]])) {
			tmpScoreList[[sampling[1]]] <- c(tmpScoreList[[sampling[1]]], sampling[2])
		} else {
			next
		}

		nuScores <- nuScore(mat, tmpScoreList)

		if (nuScores > scores && PC[sampling[1], sampling[2]] == 1) {
			scoreList <- tmpScoreList
			noChange <- 0
			scores <- nuScores
		} else {
			tmpScoreList <- scoreList
			noChange <- noChange + 1
		}
	}

	return (scoreList)

}

Score <- function(mat, PC) {
	scores <- c()
	n <- getDim(mat)
	adjMat <- matrix(0, n, n)
	nijk <- 0
	nij <- 0

	for (i in 1:n) {
		iSum <- 0
		pc <- PC[[i]]
		if (length(pc) == 0)
			next
		for (X in pc) {
			q <- getQ(mat, X)
			r <- getR(mat[, i])
			nhy <- getNhy(mat)
			wi <- getW(mat, X)
			x <- nhy/q
			y <- nhy/(r*q)
			jSum <- 0
			for (j in 1:q) {
				nij <- getNij(mat, X, wi[j], i, j, r)
				kSum <- 0
				for (k in 1:r) {
					nijk <- getNijk(mat, X, wi[j], i, j, k)
					# kSum <- kSum + log ( gamma ( nijk + y ) / gamma ( y ) )
					kSum <- kSum + lgamma ( nijk + y ) - lgamma ( y )
				}
				# jSum <- jSum + log ( gamma ( x ) / gamma ( nij + x ) ) + kSum
				jSum <- jSum + lgamma ( x ) - lgamma ( nij + x ) + kSum
			}
			# adjMat[i, X] <- jSum
			scores[i] <- jSum
			# iSum <- iSum + jSum
		}
	}

	for (i in 1:n) {
		for (j in i:n) {
			if (adjMat[i, j] == 0)
				next
			if (adjMat[i, j] >= adjMat[j, i]) {
				adjMat[i, j] <- 1
				adjMat[j, i] <- 0
			} else {
				adjMat[j, i] <- 1
				adjMat[i, j] <- 0
			}
		}
	}

	adjMat <- graph.adjacency(adjMat)

	return (adjMat)


	# n <- getDim(mat)
	# adjMat <- matrix(0, n, n)
	# nijk <- 0
	# nij <- 0
	# iSum <- 0
	# for (i in 1:n) {
	# 	pc <- PC[[i]]
	# 	if (length(pc) == 0)
	# 		next
	# 	q <- getQ(mat, pc)
	# 	r <- getR(mat[, i])
	# 	nhy <- getNhy(mat)
	# 	wi <- getW(mat, pc)
	# 	x <- nhy/q
	# 	y <- nhy/(r*q)
	# 	# print(c(q, r, nhy, x, y))
	# 	jSum <- 0
	# 	for (j in 1:q) {
	# 		if (length(pc) == 1) {
	# 			nij <- getNij(mat, pc, wi[j], i, j)
	# 		} else {
	# 			nij <- getNij(mat, pc, wi[j,], i, j)
	# 		}
	# 		kSum <- 0
	# 		for (k in 1:r) {
	# 			if (length(pc) == 1) {
	# 				nijk <- getNijk(mat, pc, wi[j], i, j, k)
	# 			} else {
	# 				nijk <- getNijk(mat, pc, wi[j,], i, j, k)
	# 			}
	# 			# print(nijk)
	# 			kSum <- kSum + log ( gamma ( (nijk/10) + y ) / gamma ( y ) )
	# 		}
	# 		jSum <- jSum + log ( gamma ( x ) / gamma ( (nij/10) + x ) ) + kSum
	# 	}
	# 	iSum <- iSum + jSum
	# }
}

# dimAdjMat <- dim(Matrixy)[2]
# AdjMatrix <<- matrix(0, 5, 5)

# aj <- graph.adjacency(mat)
# plot(aj)
# mat

# bench <- benchmark(MMPC(Matrix), mmpc(tmp), replications=1, columns = c("test", "elapsed", "relative"))

# bench <- benchmark(MaxMinHeuristic(1, 4, Matrix, c(2,5)), ForwardPhase(1, Matrix), BackwardPhase(1, 4), columns = c(1,2,3), replications = 5)
# vec <- Matrix[1,c(1,2,3)]
# cnt <- 0
# for (i in 1:dim(Matrix[,c(1,2,3)])[1]) {
# 	if(identical(vec, Matrix[i,c(1,2,3)]))
# 		cnt <- cnt + 1
# }

Test <- function() {
	# data(learning.test)
	dev.new()
	plot(mmpc(tmp))
	dev.new()
	plot(mmhc(tmp))
	myMat <- data.matrix(tmp)
	PC <- MMPC(myMat)
	dev.new()
	plot(PC)
	dev.new()
	plot(Score(myMat, PC))
}