
## Up until R 2.15.0, the require("methods") is needed but (now)
## triggers an warning from R CMD check
#.onLoad <- function(libname, pkgname){
#    #require("methods")  ## needed with R <= 2.15.0
#    loadRcppModules()
#}


## For R 2.15.1 and later this also works. Note that calling loadModule() triggers
## a load action, so this does not have to be placed in .onLoad() or evalqOnLoad().
# loadModule("Rcpp", TRUE)
# loadModule("igraph", TRUE)
loadModule("mmhc", TRUE)
source("R/data.R")

# the mmhc function
mmhc <- function(data) {
    library(Rcpp) # load Rcpp package
    library(igraph) # load igraph package
    columnNames <- colnames(df) # get the column names of the data frame
    C <- new(MMHC, data) # initalize the class object
    C$mmpc() # first reconstruct the skeleton (max-min parents and children algorithm)
    C$mmhc() # set the edges (BDeu score)
    plotObj <- C$adjMat()
    colnames(plotObj) <- columnNames
    plotObj <- graph.adjacency(plotObj) # makes a plotable object with the igraph package
    plot(plotObj) # plots the object
}