source("newCode.R")

U <- t(unique(t(Matrix)))

mmpcBench <- benchmark(MMPC(Matrix), mmpc(Example(250,char=FALSE)), replications=1, columns = c("test", "elapsed", "relative"))

StatsBench <- benchmark(Statistics(Matrix, U, card) ,replications=10, columns = c("test", "elapsed", "relative"))

solution <- MMPC(Matrix)