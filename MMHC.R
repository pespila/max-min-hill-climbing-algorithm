# library("inline")
# library("bnlearn")
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

MMPC.MightyOfSet <- function(observedMatrix) {
    n <- dim(observedMatrix)[2]
    Df <- c()
    for (i in 1:n) {
        Df[i] <- length(unique(observedMatrix[,i]))
    }
    return (Df)
}

MMPC.DegreesOfFreedom <- function(Df,target,selected,subset) {
    out <- (Df[target]-1)*(Df[selected]-1)
    if(length(subset)==1) {
        out <- out*Df[subset]
    } else {
        for (i in 1:length(subset)) {
            out <- out*Df[subset[i]]
        }
    }
    return (out)
}

MMPC.MyFunction <- function(observedMatrix,target,selected,subset) {
    out <- 0
    ABCvec <- observedMatrix[1,c(target,selected,subset)]
    ACvec <- ABCvec[c(1,3:length(ABCvec))]
    BCvec <- ABCvec[2:length(ABCvec)]
    Cvec <- ABCvec[3:length(ABCvec)]
    tmp <- 1
    partialMatrix <- observedMatrix[,c(target,selected,subset)]
    while(tmp>0) {
        ABClist <- Sabc(partialMatrix,ABCvec)
        toBeDeletedItems <- ABClist[2:length(ABClist)]
        ABCvalue <- ABClist[1]
        ACvalue <- Sabc(observedMatrix[,c(target,subset)],ACvec)[1]
        BCvalue <- Sabc(observedMatrix[,c(selected,subset)],BCvec)[1]
        Cvalue <- Sabc(as.matrix(observedMatrix[,subset]),Cvec)[1]
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
    return (out)
}

MMPC.Pvalue <- function(a,b,c) {
    tmp <- MMPC.Example(1000)
    statistic <- MMPC.MyFunction(tmp,a,b,c)
    df <- MMPC.MightyOfSet(tmp)
    df <- MMPC.DegreesOfFreedom(df,a,b,c)
    out <- pchisq(statistic, df, lower.tail = FALSE)
    return (out)
}


# MMPC.ABCvalue <- function(partialMatrix,a,b,c) {
#     partial_abc <- partialMatrix[,c(a,b,c)]
#     if(dim(partialMatrix)[1]==1) vec_abc <- partial_abc
#     else vec_abc <- partial_abc[1,]
#     out <- Sabc(vec_abc,partial_abc)
#     return (out)
# }

# MMPC.Svalue <- function(observedMatrix,a,b,c) {
#     out <- c()
#     partial_ac <- observedMatrix[,c(a,c)]
#     partial_bc <- observedMatrix[,c(b,c)]
#     if(length(c)==1) partial_c <- as.matrix(observedMatrix[,c(c)])
#     else partial_c <- observedMatrix[,c(c)]
#     vec_ac <- partial_ac[1,]
#     vec_bc <- partial_bc[1,]
#     if(length(c)==1) vec_c <- partial_c
#     else vec_c <- partial_c[1,]
#     out <- c(out, Sabc(vec_ac,partial_ac))
#     out <- c(out, Sabc(vec_bc,partial_bc))
#     out <- c(out, Sabc(vec_c,partial_c))
#     return (out)
# }

# MMPC.Gvalue <- function(observedMatrix,target,selected,subset) {
#     out <- 0
#     tmp <- 1
#     ABCval <- MMPC.ABCvalue(observedMatrix,target,selected,subset)
#     toBeDeletedItems <- c()
#     # while(tmp>0) {
#     for(i in 1:10) {
#         Svalues <- MMPC.Svalue(observedMatrix,target,selected,subset)
#         out <- out + 2*ABCval[1]*log((ABCval[1]*Svalues[3])/(Svalues[1]*Svalues[2]))
#         print(ABCval[2:length(ABCval)])
#         toBeDeletedItems <- c(toBeDeletedItems,ABCval[2:length(ABCval)])
#         tmp <- dim(observedMatrix[toBeDeletedItems,])[1]
#         ABCval <- MMPC.ABCvalue(observedMatrix[toBeDeletedItems,],target,selected,subset)
#     }
#     return (out);
# }


# MMPC.Test <- function(rows) {
#         n <- 1
#         MyMatrix <- matrix(,rows,200)
#         for(i in 1:200) {
#                 for(j in 1:rows) {
#                         MyMatrix[j,n] <- sample(c(0,1),1,replace=TRUE,prob=c(0.2,0.8))
#                 }
#                 n <- n + 1
#         }
#         rownames(MyMatrix) <- 1:rows
#         # colnames(MyMatrix) <- paste("",letters,sep="")
#         return (MyMatrix)
# }

# MMPC.Gsquared <- function(observedMatrix,target,selected,subset) {
#     # if(!length(subset)) subset <- sample(c(1:(target-1)),1)

#     # tmp <- c("1")
#     sum=0.0
#     partialMatrix <- observedMatrix[,c(target,selected,subset)]
#     # while(dim(partialMatrix)[1]>0) {
#     for(i in 1:10) {
#         rownameList <- rownames(observedMatrix[which(duplicated(observedMatrix[1,])),])
#         Sabc <- MMPC.Count(observedMatrix[,c(target,selected,subset)],observedMatrix[,c(target,selected,subset)][1,])
#         Sac <- MMPC.Count(observedMatrix[,c(target,subset)],observedMatrix[,c(target,subset)][1,])
#         Sbc <- MMPC.Count(observedMatrix[,c(selected,subset)],observedMatrix[,c(selected,subset)][1,])
#         if(length(subset)==1) Sc <- MMPC.Count(observedMatrix[,subset],observedMatrix[,subset][1],vec=TRUE)
#         else Sc <- MMPC.Count(observedMatrix[,subset],observedMatrix[,subset][1,])
#         # tmp <- c(tmp,Sabc[2:length(Sabc)])
#         print(Sabc)
#         observedMatrix <- MMPC.UpdateMatrix(observedMatrix,rownameList)
#         # print(dim(partialMatrix))
#         sum=sum+(2*Sabc*log((Sabc*Sc)/(Sac*Sbc)))
#         # print(dim(partialMatrix)[1])
#         print(sum)
#     }
#     return (sum)
# }

# MMPC.UpdateMatrix <- function(partialMatrix,toBeDeletedItems) {
#     # print(partialMatrix[!rownames(partialMatrix) %in% toBeDeletedItems[2:length(toBeDeletedItems)], ])
#     # print(toBeDeletedItems)
#     # print(dim(partialMatrix[!rownames(partialMatrix) %in% toBeDeletedItems,]))
#     return (partialMatrix[!rownames(partialMatrix) %in% toBeDeletedItems,])
# }

# MMPC.Count <- function(partialMatrix,target,vec=FALSE) {
#     # print(length(partialMatrix[which(duplicated(targetVector))])+1)
#     if(vec) tmp <- length(partialMatrix[which(duplicated(partialMatrix[1]))])+1
#     # if(vec) tmp <- c(length(partialMatrix[which(duplicated(targetVector))])+1,rownames(partialMatrix[which(duplicated(targetVector)),]))
#     else tmp <- dim(partialMatrix[which(duplicated(partialMatrix[1,])),])[1]+1
#     # else tmp <- c(dim(partialMatrix[which(duplicated(targetVector)),])[1]+1,rownames(partialMatrix[which(duplicated(targetVector)),]))
#     return (tmp)
# }



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