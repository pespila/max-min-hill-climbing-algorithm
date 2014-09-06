min <- 20
max <- 30
x<-seq(0,max,by=.25)
y<-dchisq(x,df=10)
plot(x,y,type="l")
for (i in 1:4) {
	xq=c(x[x>=min],max,min)
	yq=c(y[x>=min],0,0)
	polygon(xq,yq,col="orange")
	file <- paste('chi',i,'.eps',sep="")
	dev.copy(postscript,file)
	dev.off()
	min <- min - 5
}