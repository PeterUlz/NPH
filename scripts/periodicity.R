

#define sliding window smoothing function
slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}


findPeaks <-function (x, thresh = 0) 
{
    pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
    if (!missing(thresh)) {
        pks[x[pks - 1] - x[pks] > thresh]
    }
    else pks
}

data<-read.csv("MergedControls_Plasma_Top1000_tss.txt",header=TRUE,sep="\t")
smoothed<-slideFunct(data$Mean.Cov,50)
png("TSS_coverage_smoothed")
plot(data$Position,data$Mean.Cov,ylim=c(0,1),type="l",col="red",lwd=2)
lines(data$Position[25:1975],smoothed,col="blue")
dev.off()
findPeaks(smoothed)-975
png("Periodogram.png")
spec.pgram(fft(data$Mean.Cov[1000:2000]),spans=c(2,10))
dev.off()
pgram<-spec.pgram(fft(data$Mean.Cov[1000:2000]),spans=c(2,10),plot=FALSE)
findPeaks(pgram$spec)*1000
