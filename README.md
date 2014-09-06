max-min-hill-climbing-algorithm
===============================

This algorithm reconstructs Bayesian Networks from observational data. Therefor it first builds the skeleton of the DAG (directed acyclic graph) with the max-min parents and children (mmpc) algorithm. Afterwards it directs the edges between the vertices with the Bayesian Dirichlet likelihood-equivalence uniform score. For more information on that read the report appended or "The max-min hill-climbing Bayesian network structure learning algorithm", by Ioannis Tsamardinos, Laura E. Brown & Constantin F. Aliferis.

===============================
INSTALL
===============================
Before you can use this package, be sure you have the latest R version (>=3.0), "RCPP" version (>=0.11.1) and the "igraph" package installed.

Download the R-source file and install it into your R environment with:
install.packages("mmhc...")

Here is the example from the man pages of the package:

## Manuel Workflow: With this you are able to reconstruct the graph step by step
library(Rcpp) # load Rcpp package
library(igraph) # load igraph package
data <- student(1000) # initalize the underlying example with 1000 observations
C <- new(MMHC, data) # initalize the class object
C$mmpc() # first reconstruct the skeleton (max-min parents and children algorithm)
C$pc() # returns the PC set. It is a list where the n-th list element stands for the n-th node in your graph. The elements of one node are the parents/children of the node.
C$mmhc() # set the edges (BDeu score)
C$adjMat() # returns the adjacency matrix
C$score() # returns the score of the graph
plotObj <- graph.adjacency(C$adjMat()) # makes a plotable object with the igraph package
plot(plotObj) # plots the object

## Call the function
data <- student(1000) # as above
mmhc(data) # gives you the plot of the graph (no return value)

\name{max-min-hill-climbing-algorithm-package}
\alias{mmhc}
\docType{package}
\title{
    The max-min-hill-climbing algorithm
}
\description{
    This algorithm reconstructs Bayesian Networks from observational data. Therefor it first builds the skeleton of the DAG (directed acyclic graph) with the max-min parents and children (mmpc) algorithm. Afterwards it directs the edges between the vertices with the Bayesian Dirichlet likelihood-equivalence uniform (BDeu) score. For more information on that read the report appended or "The max-min hill-climbing Bayesian network structure learning algorithm", by Ioannis Tsamardinos, Laura E. Brown & Constantin F. Aliferis.
}
\details{
    \tabular{ll}{
    Package: \tab mmhc\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2014-09-06\cr
    License: \tab No License\cr
    
    R functions\tab ----------------------\cr
    
    student(int x) \tab returns a data frame with x observations; one of the two self made examples\cr
    rainy(int x) \tab returns a data frame with x observations; one of the two self made examples\cr
    mmhc(data.frame x) \tab input: the observed data x (a data frame); calculates everything and makes a plot\cr
    
    C++ methods \tab ----------------------\cr
    
    C <- new(MMHC, data.frame) \tab initalizes the mmhc class\cr
    C$mmpc() \tab executes the MMPC algorithm\cr
    C$mmhc() \tab executes the MMHC algorithm\cr
    C$pc() \tab this C++-methods returns the PC set\cr
    C$adjMat() \tab this C++-methods returns the adjacency matrix\cr
    C$score() \tab this C++-methods returns the score of the reconstructed graph\cr
    C$mat() \tab this C++-methods returns the data frame converted into an integer matrix\cr
    }
}
\author{
    Maintainer: Michael Bauer <michael1.bauer@ur.de>
}
\references{
    "The max-min hill-climbing Bayesian network structure learning algorithm", by Ioannis Tsamardinos, Laura E. Brown & Constantin F. Aliferis.\cr
    "A Scoring Function for Learning Bayesian Networks based on Mutual Information and
    Conditional Independence Tests", by Luis M. de Campos.
}
\keyword{ package }
\seealso{
    http://www.dsl-lab.org/supplements/mmhc_paper/paper_online.pdf\cr
    bnlearn-package
}
\examples{
    ## Call the function
    
    data <- student(1000) # as above
    mmhc(data) # gives you the plot of the graph (no return value)
    
    
    ## Manuel Workflow: With this you are able to reconstruct the graph step by step
    
    library(Rcpp) # load Rcpp package
    library(igraph) # load igraph package
    data <- student(1000) # initalize the underlying example with 1000 observations
    C <- new(MMHC, data) # initalize the class object
    C$mmpc() # first reconstruct the skeleton (max-min parents and children algorithm)
    C$pc() # returns the PC set. It is a list where the n-th list element stands for the n-th node in your graph. The elements of one node are the parents/children of the node.
    C$mmhc() # set the edges (BDeu score)
    C$adjMat() # returns the adjacency matrix
    C$score() # returns the score of the graph
    plotObj <- graph.adjacency(C$adjMat()) # makes a plotable object with the igraph package
    plot(plotObj) # plots the object
}
