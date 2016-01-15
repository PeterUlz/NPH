data<-read.csv("MergedControls_allGenes_prediction.csv",header=TRUE,sep="\t")

#Check single parameters
png("Plots/Broad_Coverage_hist.png")
hist(data$Broad.TSS.Coverage,col=rgb(1,0,0,0.5),xlim=c(0,30),breaks=100)
dev.off()
png("Plots/Broad_Coverage_normalized_hist.png")
hist(data$Broad.TSS.Coverage.Norm,col=rgb(1,0,0,0.5),xlim=c(0,1.5),breaks=50)
dev.off()
png("Plots/Small_Coverage_hist.png")
hist(data$Small.TSS.coverage,col=rgb(0.2,0.8,0.3,0.5),xlim=c(0,2),breaks=100)
dev.off()
png("Plots/5Prime_Slope_hist.png")
hist(data$X5.Slope,freq=FALSE,col=rgb(1,0,0,0.5),xlim=c(-0.03,0.03),breaks=150)
dev.off()
png("Plots/3Prime_Slope_hist.png")
hist(data$X3.Slope,freq=FALSE,col=rgb(1,0,0,0.5),xlim=c(-0.03,0.03),breaks=50)
dev.off()

#Additional Plots
png("Plots/Broad_vs_Small_TSS_Coverage.png")
plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
dev.off()
png("Plots/5prime_vs_3prime_slope.png")
plot(data$X5.Slope,data$X3.Slope,xlab="5' Slope",ylab="3' Slope",xlim=c(-0.02,0.02),col=cols,ylim=c(-0.02,0.02))
dev.off()

#############################################################################################################################################################################
# Analysis functions to determine clustering versions
#
#
#
#############################################################################################################################################################################
analyze_and_plot_clustering<-function(data,output_dir) {
    cols<-rep("black",length(data$Type))
    cols[which(data$Clustering == expressed_cluster)]<-"red"
    cols[which(data$Clustering != expressed_cluster)]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_clustering.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()

    y<-which(data$Clustering != expressed_cluster)
    gene_names<-data.frame(data$Gene[y], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

    x<-which(data$Clustering == expressed_cluster)
    gene_names<-data.frame(data$Gene[x], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

}
#############################################################################################################################################################################
sensitivity_top1000<-function(data,output_dir){
    #Check in more detail whether Top1000 genes are in expressde cluster and bottom 1000 in unexpressed cluster:
    top1000<-read.csv("../../../ref/Plasma-RNASeq/Top1000_NMonly.txt",header=FALSE)
    bottom1000<-read.csv("../../../ref/Plasma-RNASeq/Bottom1000_NMonly.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    data$Expressed<-rep("unknown",length(data$Gene))
    for (top in top1000$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        data$Expressed[index] <- "Top1000"
    }
    for (bottom in bottom1000$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        data$Expressed[index] <- "Bottom1000"
    }
    true_expressed<-which(data$Expressed == "Top1000" & data$Clustering == expressed_cluster)
    true_unexpressed<-which(data$Expressed == "Bottom1000" & data$Clustering != expressed_cluster)
    false_expressed<-which(data$Expressed == "Bottom1000" & data$Clustering == expressed_cluster)
    false_unexpressed<-which(data$Expressed == "Top1000" & data$Clustering != expressed_cluster)

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(data$Expressed == "Top1000" & data$Clustering != expressed_cluster)]<-"red"
    cols[which(data$Expressed == "Bottom1000" & data$Clustering == expressed_cluster)]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    points(data$Broad.TSS.Coverage.Norm[true_expressed],data$Small.TSS.coverage[true_expressed],,col=cols[true_expressed])
    points(data$Broad.TSS.Coverage.Norm[true_unexpressed],data$Small.TSS.coverage[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad.TSS.Coverage.Norm[false_unexpressed],data$Small.TSS.coverage[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad.TSS.Coverage.Norm[false_expressed],data$Small.TSS.coverage[false_expressed],,col=cols[false_expressed])
    dev.off()

    gene_names<-data.frame(data$Gene[true_expressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Top1000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[true_unexpressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Bottom1000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[false_unexpressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Top1000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[false_expressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Bottom1000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_housekeeping<-function(data,output_dir){
    #Check in more detail whether Housekeeping genes are assigned into expressed cluster:
    housekeeping<-read.csv("../../../ref/Housekeeping/HK_gene_names.txt",header=FALSE)

    housekeeping_in_expressed_cluster<-0
    housekeeping_in_unexpressed_cluster<-0
    data$Housekeeping<-rep("unknown",length(data$Gene))
    for (top in housekeeping$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        data$Housekeeping[index] <- "Housekeeping"
    }
    true_expressed<-which(data$Housekeeping == "Housekeeping" & data$Clustering == expressed_cluster)
    false_unexpressed<-which(data$Housekeeping == "Housekeeping" & data$Clustering != expressed_cluster)

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Housekeeping.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()

    gene_names<-data.frame(data$Gene[true_expressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_Housekeeping_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[false_unexpressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_Housekeeping_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")


    write.table(t(c("Housekeeping (Eisenberg) in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Houskeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Housekeeping (Eisenberg) in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Houskeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}

#############################################################################################################################################################################
clustering_vs_rma<-function(data,output_dir){
    #check RMA values for genes and plot "Heatmap"
    data$rma<-rep(0,length(data$Gene))
    rma_values<-read.csv("../../../ref/Plasma-RNASeq/NonPregnant_annotated_noChrM_dedup_onlyNM.txt.csv",header=TRUE,sep="\t")
    for (gene in rma_values$Gene) {
        index<-which(data$Gene == gene)
        rma_index<-which(rma_values$Gene == gene)
        if (length(index) == 0) {
            next
        }
        if (length(rma_index) == 1) {
            rma_value <- rma_values$Mean[rma_index]
        }
        else  {
            rma_value <- mean(rma_values$Mean[rma_index])
        }
        data$rma[index] <- rep(rma_value,length(index))
    }
    colors<-colorRampPalette(c("green","red"))(3)
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    level<-which(data$rma > 7)
    cols[level]<-colors[3]
    level<-which(data$rma < 4)
    cols[level]<-colors[1]
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Heatmap.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()

    #compare distance to expressed cluster center to RMA
    colors<-colorRampPalette(c("green","red"))(max(floor(data$rma)))
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    for (i in min(floor(data$rma)):max(floor(data$rma))) {
       level<-which(floor(data$rma) == i)
       cols[level]<-colors[i]
    }
    png(paste(output_dir,"/Distance_vs_RMA.png",sep=""))
    data$distance_expressed<-sqrt((data$Broad.TSS.Coverage.Norm - clustering$centers[expressed_cluster,3])^2 + (data$Small.TSS.coverage - clustering$centers[expressed_cluster,4])^2 + (data$X5.Slope - clustering$centers[expressed_cluster,1])^2 +(data$X3.Slope - clustering$centers[expressed_cluster,2])^2 )
    plot(data$distance_expressed,data$rma,col=cols,xlim=c(0,1))
    dev.off()
}

#############################################################################################################################################################################
clustering_vs_fpkm<-function(data,output_dir) {
    #check FPKM values for genes and plot "Heatmap"
    data$fpkm<-rep(0,length(data$Gene))
    fpkm_values<-read.csv("../../../ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",header=TRUE,sep="\t")
    for (gene in fpkm_values$tracking_id) {
        index<-which(data$Gene == gene)
        fpkm_index<-which(fpkm_values$tracking_id == gene)
        if (length(index) == 0) {
            next
        }
        if (length(fpkm_index) == 1) {
            fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
        }
        else  {
            fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
        }
        data$fpkm[index] <- rep(fpkm_value,length(index))
    }
    true_expressed<-which(data$fpkm > 8 & data$Clustering == expressed_cluster)
    true_unexpressed<-which(data$fpkm > 0 & data$fpkm < 0.1 & data$Clustering != expressed_cluster)
    false_expressed<-which(data$fpkm > 0 & data$fpkm < 0.1 & data$Clustering ==  expressed_cluster)
    false_unexpressed<-which(data$fpkm > 8 & data$Clustering != expressed_cluster)

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(0.4,0,0.77)
    cols[true_unexpressed]<-rgb(0.4,0.6,1)
    cols[false_unexpressed]<-rgb(1,0,0)
    cols[false_expressed]<-rgb(1,0,0)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Heatmap.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Highlight.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    points(data$Broad.TSS.Coverage.Norm[false_unexpressed],data$Small.TSS.coverage[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad.TSS.Coverage.Norm[false_expressed],data$Small.TSS.coverage[false_expressed],,col=cols[false_expressed])
    points(data$Broad.TSS.Coverage.Norm[true_expressed],data$Small.TSS.coverage[true_expressed],,col=cols[true_expressed])
    points(data$Broad.TSS.Coverage.Norm[true_unexpressed],data$Small.TSS.coverage[true_unexpressed],,col=cols[true_unexpressed])
    dev.off()


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_top100<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv("../../../ref/Plasma-RNASeq/Top100_NMonly.txt",header=FALSE)
    bottom100<-read.csv("../../../ref/Plasma-RNASeq/Bottom100_NMonly.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    data$Expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        data$Expressed[index] <- "Top100"
    }
    for (bottom in bottom100$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        data$Expressed[index] <- "Bottom100"
    }
    true_expressed<-which(data$Expressed == "Top100" & data$Clustering == expressed_cluster)
    true_unexpressed<-which(data$Expressed == "Bottom100" & data$Clustering != expressed_cluster)
    false_expressed<-which(data$Expressed == "Bottom100" & data$Clustering == expressed_cluster)
    false_unexpressed<-which(data$Expressed == "Top100" & data$Clustering != expressed_cluster)

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(data$Expressed == "Top100" & data$Clustering != expressed_cluster)]<-"red"
    cols[which(data$Expressed == "Bottom100" & data$Clustering == expressed_cluster)]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top100.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100.png",sep=""))
    plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
    points(data$Broad.TSS.Coverage.Norm[true_expressed],data$Small.TSS.coverage[true_expressed],,col=cols[true_expressed])
    points(data$Broad.TSS.Coverage.Norm[true_unexpressed],data$Small.TSS.coverage[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad.TSS.Coverage.Norm[false_unexpressed],data$Small.TSS.coverage[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad.TSS.Coverage.Norm[false_expressed],data$Small.TSS.coverage[false_expressed],,col=cols[false_expressed])
    dev.off()


    gene_names<-data.frame(data$Gene[true_expressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[true_unexpressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[false_unexpressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

    gene_names<-data.frame(data$Gene[false_expressed], stringsAsFactors=FALSE)
    names(gene_names)<-c("GeneName")
    write.table(gene_names$GeneName, file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)

}

##### Clustering #####################################################################################
#1) First try clustering with four parameters K-means

#K-means clustering using 4 parameters
data_numeric<-rbind(data$X5.Slope,data$X3.Slope,data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage)
clustering<-kmeans(t(data_numeric),2)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1,3] > clustering$centers[2,3]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"1_Kmeans_4param")
clustering_vs_fpkm(data,"1_Kmeans_4param")
clustering_vs_rma(data,"1_Kmeans_4param")
sensitivity_top1000(data,"1_Kmeans_4param")
sensitivity_top100(data,"1_Kmeans_4param")
##########################################################################################
#2) Clustering with two parameters K-means (Broad and Small TSS coverage)
#K-means clustering using 2 parameters
data_numeric<-rbind(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage)
clustering<-kmeans(t(data_numeric),2)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1,1] > clustering$centers[2,1]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"2_Kmeans_2param")
clustering_vs_fpkm(data,"2_Kmeans_2param")
clustering_vs_rma(data,"2_Kmeans_2param")
sensitivity_top1000(data,"2_Kmeans_2param")
sensitivity_top100(data,"2_Kmeans_2param")
##########################################################################################
#3) Clustering with four parameters K-means, but scale parameters

#K-means clustering using 4 parameters
data_numeric<-rbind(data$X5.Slope,data$X3.Slope,data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage)
data_scaled<-scale(data_numeric)
clustering<-kmeans(t(data_scaled),2)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1,3] > clustering$centers[2,3]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"3_Kmeans_4param_scaled")
clustering_vs_fpkm(data,"3_Kmeans_4param_scaled")
clustering_vs_rma(data,"3_Kmeans_4param_scaled")
sensitivity_top1000(data,"3_Kmeans_4param_scaled")
sensitivity_top100(data,"3_Kmeans_4param_scaled")
##########################################################################################
#4) Clustering with two parameters K-means, use a priori cluster centers

#K-means clustering using 2 parameters
data_numeric<-rbind(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage)
expressed<-c(0.5,0.5)
unexpressed<-c(1,1)
starting_centers<-rbind(expressed,unexpressed)
clustering<-kmeans(t(data_numeric),starting_centers)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1,1] > clustering$centers[2,1]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"4_Kmeans_explicit_centers")
clustering_vs_fpkm(data,"4_Kmeans_explicit_centers")
clustering_vs_rma(data,"4_Kmeans_explicit_centers")
sensitivity_top1000(data,"4_Kmeans_explicit_centers")
sensitivity_housekeeping(data,"4_Kmeans_explicit_centers")
sensitivity_top100(data,"4_Kmeans_explicit_centers")
######################################################################################
#5) Clustering with only one parameter K-means (Broad TSS coverage)

#K-means clustering using 1 parameters
data_numeric<-rbind(data$Broad.TSS.Coverage.Norm)
clustering<-kmeans(t(data_numeric),2)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1] > clustering$centers[2]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"5_Kmeans_onlyBroadTSSCoverage")
clustering_vs_fpkm(data,"5_Kmeans_onlyBroadTSSCoverage")
clustering_vs_rma(data,"5_Kmeans_onlyBroadTSSCoverage")
sensitivity_top1000(data,"5_Kmeans_onlyBroadTSSCoverage")
sensitivity_housekeeping(data,"5_Kmeans_onlyBroadTSSCoverage")
sensitivity_top100(data,"5_Kmeans_onlyBroadTSSCoverage")
##########################################################################################
#5.5) Clustering with two parameters K-means, use a priori cluster centers, weight Broad Coverage twice as much

#K-means clustering using 2 parameters
data_numeric<-rbind((data$Broad.TSS.Coverage.Norm*2),data$Small.TSS.coverage)
expressed<-c(0.5,0.5)
unexpressed<-c(1,1)
starting_centers<-rbind(expressed,unexpressed)
clustering<-kmeans(t(data_numeric),starting_centers)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1,1] > clustering$centers[2,1]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"5_5_Kmeans_2xBroad_1xSmall")
clustering_vs_fpkm(data,"5_5_Kmeans_2xBroad_1xSmall")
clustering_vs_rma(data,"5_5_Kmeans_2xBroad_1xSmall")
sensitivity_top1000(data,"5_5_Kmeans_2xBroad_1xSmall")
sensitivity_housekeeping(data,"5_5_Kmeans_2xBroad_1xSmall")
sensitivity_top100(data,"5_5_Kmeans_2xBroad_1xSmall")
##########################################################################################
#6) Clustering with only one parameter K-means (Small TSS coverage)

#K-means clustering using 1 parameters
data_numeric<-rbind(data$Small.TSS.Coverage)
clustering<-kmeans(t(data_numeric),2)
data$Clustering<-clustering$cluster

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1] > clustering$centers[2]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"6_Kmeans_onlySmallTSSCoverage")
clustering_vs_fpkm(data,"6_Kmeans_onlySmallTSSCoverage")
clustering_vs_rma(data,"6_Kmeans_onlySmallTSSCoverage")
sensitivity_top1000(data,"6_Kmeans_onlySmallTSSCoverage")
sensitivity_top100(data,"6_Kmeans_onlySmallTSSCoverage")

##########################################################################################
#7) Clustering with two parameter PAM (partitioning around medoids) (Broad and Small TSS coverage)
library(cluster)
#K-means clustering using 2 parameters
data_numeric<-rbind(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage)
clustering<-pam(t(data_numeric),2)
data$Clustering<-clustering$clustering

#search for expressed cluster (should have lower mean coverage)
expressed_cluster<-1
if (clustering$centers[1,1] > clustering$centers[2,1]) {
    expressed_cluster<-2
}

analyze_and_plot_clustering(data,"7_PAM_2param")
clustering_vs_fpkm(data,"7_PAM_2param")
clustering_vs_rma(data,"7_PAM_2param")
sensitivity_top1000(data,"7_PAM_2param")
sensitivity_top100(data,"7_PAM_2param")

