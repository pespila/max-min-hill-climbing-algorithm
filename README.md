max-min-hill-climbing-algorithm
===============================

This algorithm reconstructs Bayesian Networks from observational data. Therefor it first builds the skeleton of the DAG (directed acyclic graph) with the max-min parents and children (mmpc) algorithm. Afterwards it directs the edges between the vertices with the Bayesian Dirichlet likelihood-equivalence uniform score. For more information on that read the report appended or "The max-min hill-climbing Bayesian network structure learning algorithm", by Ioannis Tsamardinos, Laura E. Brown & Constantin F. Aliferis.

INSTALL
-------------------------------
Before you can use this package, be sure you have the latest R version (>=3.0), "RCPP" version (>=0.11.1) and the "igraph" package installed.

Download the R-source file and install it into your R environment with:

install.packages("mmhc...")

Here is the example from the man pages of the package:

## Basic R functions
    
+ **student(int x)**   *returns a data frame with x observations; one of the two self made examples*
+ **rainy(int x)**   *returns a data frame with x observations; one of the two self made examples*
+ **mmhc(data.frame x)**   *input: the observed data x (a data frame); calculates everything and makes a plot*
    
## Basic C++ methods
    
+ **C <- new(MMHC, data.frame)**   *initalizes the mmhc class*
+ **C$mmpc()**   *executes the MMPC algorithm*
+ **C$mmhc()**   *executes the MMHC algorithm*
+ **C$pc()**   *this C++-methods returns the PC set*
+ **C$adjMat()**   *this C++-methods returns the adjacency matrix*
+ **C$score()**   *this C++-methods returns the score of the reconstructed graph*
+ **C$mat()**   *this C++-methods returns the data frame converted into an integer matrix*

How to use the algorithm
-------------------------------

+ **data<- student(1000)** *# as above*
+ **mmhc(data)** *# gives you the plot of the graph (no return value)*
    
Manuel Workflow
-------------------------------

### Producing the data step by step
    
+ **library(Rcpp)** *# load Rcpp package*
+ **library(igraph)** *# load igraph package*
+ **data <- student(1000)** *# initalize the underlying example with 1000 observations*
+ **C <- new(MMHC, data)** *# initalize the class object*
+ **C$mmpc()** *# first reconstruct the skeleton (max-min parents and children algorithm)*
+ **C$pc()** *# returns the PC set. It is a list where the n-th list element stands for the n-th node in your graph. The elements of one node are the parents/children of the node.*
+ **C$mmhc()** *# set the edges (BDeu score)*
+ **C$adjMat()** *# returns the adjacency matrix*
+ **C$score()** *# returns the score of the graph*
+ **plotObj <- graph.adjacency(C$adjMat())** *# makes a plotable object with the igraph package*
+ **plot(plotObj)** *# plots the object*