library("inline")
library("bnlearn")
library("Rcpp")
sourceCpp("mmpc.cpp");

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
        # observedMatrix <- data.frame(observedMatrix)
        return (observedMatrix)
}

MMPC.Test <- function(rows) {
        n <- 1
        MyMatrix <- matrix(,rows,200)
        for(i in 1:200) {
                for(j in 1:rows) {
                        MyMatrix[j,n] <- sample(c(0,1),1,replace=TRUE,prob=c(0.2,0.8))
                }
                n <- n + 1
        }
        rownames(MyMatrix) <- 1:rows
        # colnames(MyMatrix) <- paste("",letters,sep="")
        return (MyMatrix)
}

MMPC.Gsquared <- function(observedMatrix,target,selected,subset) {
    # if(!length(subset)) subset <- sample(c(1:(target-1)),1)

    # tmp <- c("1")
    sum=0.0
    partialMatrix <- observedMatrix[,c(target,selected,subset)]
    # while(dim(partialMatrix)[1]>0) {
    for(i in 1:10) {
        rownameList <- rownames(observedMatrix[which(duplicated(observedMatrix[1,])),])
        Sabc <- MMPC.Count(observedMatrix[,c(target,selected,subset)],observedMatrix[,c(target,selected,subset)][1,])
        Sac <- MMPC.Count(observedMatrix[,c(target,subset)],observedMatrix[,c(target,subset)][1,])
        Sbc <- MMPC.Count(observedMatrix[,c(selected,subset)],observedMatrix[,c(selected,subset)][1,])
        if(length(subset)==1) Sc <- MMPC.Count(observedMatrix[,subset],observedMatrix[,subset][1],vec=TRUE)
        else Sc <- MMPC.Count(observedMatrix[,subset],observedMatrix[,subset][1,])
        # tmp <- c(tmp,Sabc[2:length(Sabc)])
        print(Sabc)
        observedMatrix <- MMPC.UpdateMatrix(observedMatrix,rownameList)
        # print(dim(partialMatrix))
        sum=sum+(2*Sabc*log((Sabc*Sc)/(Sac*Sbc)))
        # print(dim(partialMatrix)[1])
        print(sum)
    }
    return (sum)
}

MMPC.UpdateMatrix <- function(partialMatrix,toBeDeletedItems) {
    # print(partialMatrix[!rownames(partialMatrix) %in% toBeDeletedItems[2:length(toBeDeletedItems)], ])
    # print(toBeDeletedItems)
    # print(dim(partialMatrix[!rownames(partialMatrix) %in% toBeDeletedItems,]))
    return (partialMatrix[!rownames(partialMatrix) %in% toBeDeletedItems,])
}

MMPC.Count <- function(partialMatrix,target,vec=FALSE) {
    # print(length(partialMatrix[which(duplicated(targetVector))])+1)
    if(vec) tmp <- length(partialMatrix[which(duplicated(partialMatrix[1]))])+1
    # if(vec) tmp <- c(length(partialMatrix[which(duplicated(targetVector))])+1,rownames(partialMatrix[which(duplicated(targetVector)),]))
    else tmp <- dim(partialMatrix[which(duplicated(partialMatrix[1,])),])[1]+1
    # else tmp <- c(dim(partialMatrix[which(duplicated(targetVector)),])[1]+1,rownames(partialMatrix[which(duplicated(targetVector)),]))
    return (tmp)
}

tmp=MMPC.Example(100)



# print(sexp_type(tmp))
# dyn.load("main.so")
# print(.Call("MyMatrix",tmp))

# sexp_type <- cfunction(c(x = "ANY"), '
#   switch (TYPEOF(x)) {
#     case NILSXP:      return mkString("NILSXP");
#     case SYMSXP:      return mkString("SYMSXP");
#     case LISTSXP:     return mkString("LISTSXP");
#     case CLOSXP:      return mkString("CLOSXP");
#     case ENVSXP:      return mkString("ENVSXP");
#     case PROMSXP:     return mkString("PROMSXP");
#     case LANGSXP:     return mkString("LANGSXP");
#     case SPECIALSXP:  return mkString("SPECIALSXP");
#     case BUILTINSXP:  return mkString("BUILTINSXP");
#     case CHARSXP:     return mkString("CHARSXP");
#     case LGLSXP:      return mkString("LGLSXP");
#     case INTSXP:      return mkString("INTSXP");
#     case REALSXP:     return mkString("REALSXP");
#     case CPLXSXP:     return mkString("CPLXSXP");
#     case STRSXP:      return mkString("STRSXP");
#     case DOTSXP:      return mkString("DOTSXP");
#     case ANYSXP:      return mkString("ANYSXP");
#     case VECSXP:      return mkString("VECSXP");
#     case EXPRSXP:     return mkString("EXPRSXP");
#     case BCODESXP:    return mkString("BCODESXP");
#     case EXTPTRSXP:   return mkString("EXTPTRSXP");
#     case WEAKREFSXP:  return mkString("WEAKREFSXP");
#     case S4SXP:       return mkString("S4SXP");
#     case RAWSXP:      return mkString("RAWSXP");
#     default:          return mkString("<unknown>");
# }')