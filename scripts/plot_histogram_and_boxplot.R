data<-read.csv(args[4],header=TRUE,sep="\t")
expr<-which(data$Type == "Expressed")
unexpr<-which(data$Type == "Unexpressed")
boxplot(data$Coverage[expr],data$Coverage[unexpr],ylim=c(0,20),col=c("red","green"))
hist(data$Coverage[unexpr],col=rgb(0.2,0.8,0.3,0.5),xlim=c(0,50),breaks=100)
hist(data$Coverage[expr],col=rgb(1,0,0,0.5),xlim=c(0,50),add=TRUE,breaks=100)

boxplot(data$X5.Slope[expr],data$X5.Slope[unexpr],ylim=c(-0.01,0.01),col=c("red","green"))
hist(data$X5.Slope[expr],freq=FALSE,col=rgb(1,0,0,0.5),xlim=c(-0.03,0.03),breaks=50)
hist(data$X5.Slope[unexpr],freq=FALSE,col=rgb(0.2,0.8,0.3,0.5),xlim=c(-0.03,0.03),breaks=2000,add=TRUE)

boxplot(data$X3.Slope[expr],data$X3.Slope[unexpr],ylim=c(-0.01,0.01),col=c("red","green"))
hist(data$X3.Slope[expr],freq=FALSE,col=rgb(1,0,0,0.5),xlim=c(-0.01,0.01),breaks=100)
hist(data$X3.Slope[unexpr],freq=FALSE,col=rgb(0.2,0.8,0.3,0.5),xlim=c(-0.01,0.01),breaks=2000,add=TRUE)

boxplot(data$Amplitude[expr],data$Amplitude[unexpr],ylim=c(-0.01,0.05),col=c("red","green"))
hist(data$Amplitude[expr],col=rgb(1,0,0,0.5),xlim=c(-0.01,0.05),breaks=100)
hist(data$Amplitude[unexpr],col=rgb(0.2,0.8,0.3,0.5),xlim=c(-0.01,0.05),breaks=100,add=TRUE)

boxplot(data$Coverage[expr],data$Coverage[unexpr],ylim=c(0,2),col=c("red","green"))
hist(data$Coverage[unexpr],col=rgb(0.2,0.8,0.3,0.5),xlim=c(0,2),breaks=20)
hist(data$Coverage[expr],col=rgb(1,0,0,0.5),xlim=c(0,2),add=TRUE,breaks=20)


data<-read.csv("output/PredictActiveGenes/TSS_coverage_dilution_series.txt",header=TRUE,sep="\t")
png("Dilution_Series.barplot.png")
barx<-barplot(data$Mean,names.arg=data$Geneset,las=2,col=colorRampPalette(c("green","red"))(9),ylim=c(0,20))
arrows(barx,data$UpperBound, barx, data$Mean, angle=90, code=1)
arrows(barx,data$LowerBound, barx, data$Mean, angle=90, code=1)
dev.off()
