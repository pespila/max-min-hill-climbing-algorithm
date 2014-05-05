if(!("Rcpp" %in% rownames(installed.packages()))) install.packages("Rcpp")
if(!("sets" %in% rownames(installed.packages()))) install.packages("sets")

require("Rcpp")
require("sets")
sourceCpp("mmhc.cpp");

MMPC <- function() UseMethod("MMPC")

MMPC.Example <- function(rows,cols=5) {
        names <- c("Difficulty","Intelligence","SAT","Grade","Letter")
        n <- 1
        observedMatrix <- matrix(,rows,cols)
        observedMatrix[,1] <- sample(c(0,1),rows,replace=TRUE,prob=c(0.6,0.4))
        observedMatrix[,2] <- sample(c(0,1),rows,replace=TRUE,prob=c(0.7,0.3))
        for (i in observedMatrix[,2]) {
                if(!i) observedMatrix[n,3] <- sample(c(0,1),1,replace=TRUE,prob=c(0.95,0.05))
                if(i) observedMatrix[n,3] <- sample(c(0,1),1,replace=TRUE,prob=c(0.2,0.8))
                if(i&&observedMatrix[n,1]) observedMatrix[n,4] <- sample(c(1,2,3),1,replace=TRUE,prob=c(0.5,0.3,0.2))
                if(i&&!observedMatrix[n,1]) observedMatrix[n,4] <- sample(c(1,2,3),1,replace=TRUE,prob=c(0.9,0.08,0.02))
                if(!i&&observedMatrix[n,1]) observedMatrix[n,4] <- sample(c(1,2,3),1,replace=TRUE,prob=c(0.05,0.25,0.7))
                if(!i&&!observedMatrix[n,1]) observedMatrix[n,4] <- sample(c(1,2,3),1,replace=TRUE,prob=c(0.3,0.4,0.3))
                if(observedMatrix[n,4]) observedMatrix[n,5] <- sample(c(0,1),1,replace=TRUE,prob=c(0.1,0.9))
                if(observedMatrix[n,4]==2) observedMatrix[n,5] <- sample(c(0,1),1,replace=TRUE,prob=c(0.4,0.6))
                if(observedMatrix[n,4]==3) observedMatrix[n,5] <- sample(c(0,1),1,replace=TRUE,prob=c(0.99,0.01))
                n <- n+1
        }
        rownames(observedMatrix) <- 1:rows
        colnames(observedMatrix) <- names
        return (observedMatrix)
}

MMPC.Example2 <- function(rows,cols=5) {
        names <- c("Difficulty","Intelligence","SAT","Grade","Letter")
        n <- 1
        observedMatrix <- matrix(,rows,cols)
        observedMatrix[,1] <- sample(c("high","low"),rows,replace=TRUE,prob=c(0.6,0.4))
        observedMatrix[,2] <- sample(c("high","low"),rows,replace=TRUE,prob=c(0.7,0.3))
        for (i in observedMatrix[,2]) {
                if(i=="low") observedMatrix[n,3] <- sample(c("bad","good"),1,replace=TRUE,prob=c(0.95,0.05))
                if(i=="high") observedMatrix[n,3] <- sample(c("bad","good"),1,replace=TRUE,prob=c(0.2,0.8))
                if(i=="high"&&observedMatrix[n,1]=="high") observedMatrix[n,4] <- sample(c("A","B","C"),1,replace=TRUE,prob=c(0.5,0.3,0.2))
                if(i=="high"&&observedMatrix[n,1]=="low") observedMatrix[n,4] <- sample(c("A","B","C"),1,replace=TRUE,prob=c(0.9,0.08,0.02))
                if(i=="low"&&observedMatrix[n,1]=="high") observedMatrix[n,4] <- sample(c("A","B","C"),1,replace=TRUE,prob=c(0.05,0.25,0.7))
                if(i=="low"&&observedMatrix[n,1]=="low") observedMatrix[n,4] <- sample(c("A","B","C"),1,replace=TRUE,prob=c(0.3,0.4,0.3))
                if(observedMatrix[n,4]=="A") observedMatrix[n,5] <- sample(c("bad","good"),1,replace=TRUE,prob=c(0.1,0.9))
                if(observedMatrix[n,4]=="B") observedMatrix[n,5] <- sample(c("bad","good"),1,replace=TRUE,prob=c(0.4,0.6))
                if(observedMatrix[n,4]=="C") observedMatrix[n,5] <- sample(c("bad","good"),1,replace=TRUE,prob=c(0.99,0.01))
                n <- n+1
        }
        rownames(observedMatrix) <- 1:rows
        colnames(observedMatrix) <- names
        return (observedMatrix)
}

MMPC.Cardinality <- function(observedMatrix) {
    n <- dim(observedMatrix)[2]
    Df <- c()
    for (i in 1:n) {
        Df[i] <- length(unique(observedMatrix[,i]))
    }
    return (Df)
}

# MMPC.DegreesOfFreedom <- function(Df,target,selected,subset) {
#     out <- (Df[target]-1)*(Df[selected]-1)
#     if(length(subset)==1) {
#         out <- out*Df[subset]
#     } else {
#         for (i in 1:length(subset)) {
#             out <- out*Df[subset[i]]
#         }
#     }
#     return (out)
# }

MMPC.Statistic <- function(observedMatrix,target,selected,subset) {
    out <- 0
    ABCvec <- observedMatrix[1,c(target,selected,subset)]
    ACvec <- ABCvec[c(1,3:length(ABCvec))]
    BCvec <- ABCvec[2:length(ABCvec)]
    Cvec <- ABCvec[3:length(ABCvec)]
    tmp <- 1
    partialMatrix <- observedMatrix[,c(target,selected,subset)]
    type <- 1
    if(class(observedMatrix[1,1])=="character") type<-0
    while(tmp>0) {
        ABClist <- Test(partialMatrix,ABCvec,type)
        toBeDeletedItems <- ABClist[2:length(ABClist)]
        ABCvalue <- ABClist[1]
        ACvalue <- Test(observedMatrix[,c(target,subset)],ACvec,type)[1]
        BCvalue <- Test(observedMatrix[,c(selected,subset)],BCvec,type)[1]
        Cvalue <- Test(as.matrix(observedMatrix[,subset]),Cvec,type)[1]
        partialMatrix <- partialMatrix[toBeDeletedItems,]
        if(length(partialMatrix)==length(c(target,selected,subset))) tmp <- 1
        else tmp <- dim(as.matrix(partialMatrix))[1]
        if(tmp>0) {
            if(length(partialMatrix)==length(c(target,selected,subset))) {
                    ABCvec <- partialMatrix
                    tmp <- 0
                }
            else ABCvec <- partialMatrix[1,]
            ACvec <- ABCvec[c(1,3:length(ABCvec))]
            BCvec <- ABCvec[2:length(ABCvec)]
            Cvec <- ABCvec[3:length(ABCvec)]
        }
        out <- out + 2*ABCvalue*log((ABCvalue*Cvalue)/(ACvalue*BCvalue))
    }
    TOL <- 0.001
    if(out < TOL) out <- 0
    return (out)
}

MMPC.Pvalue <- function(a,b,c) {
    tmp <- MMPC.Example2(1000)
    statistic <- MMPC.Statistic(tmp,a,b,c)
    df <- MMPC.Cardinality(tmp)
    df <- Df(df,a,b,c)
    out <- pchisq(statistic, df, lower.tail = FALSE)
    return (out)
}