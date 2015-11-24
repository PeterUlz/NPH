data<-read.csv(args[4],header=TRUE,sep="\t")
expr<-which(data$Type == "Expressed")
unexpr<-which(data$Type == "Unexpressed")
boxplot(data$Coverage[expr],data$Coverage[unexpr],ylim=c(0,20),col=c("red","green"))
hist(data$Coverage[unexpr],col=rgb(0.2,0.8,0.3,0.5),xlim=c(0,20),breaks=10000)
hist(data$Coverage[expr],col=rgb(1,0,0,0.5),xlim=c(0,20),add=TRUE,breaks=1000)

boxplot(data$X5.Slope[expr],data$X5.Slope[unexpr],ylim=c(-0.01,0.01),col=c("red","green"))
hist(data$X5.Slope[expr],freq=FALSE,col=rgb(1,0,0,0.5),xlim=c(-0.01,0.01),breaks=50)
hist(data$X5.Slope[unexpr],freq=FALSE,col=rgb(0.2,0.8,0.3,0.5),xlim=c(-0.01,0.01),breaks=2000,add=TRUE)

boxplot(data$X3.Slope[expr],data$X3.Slope[unexpr],ylim=c(-0.01,0.01),col=c("red","green"))
hist(data$X3.Slope[expr],freq=FALSE,col=rgb(1,0,0,0.5),xlim=c(-0.01,0.01),breaks=100)
hist(data$X3.Slope[unexpr],freq=FALSE,col=rgb(0.2,0.8,0.3,0.5),xlim=c(-0.01,0.01),breaks=2000,add=TRUE)
