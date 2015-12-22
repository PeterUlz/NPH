library(e1071)
data<-read.csv("MergedControls_allGenes_prediction.csv",header=TRUE,sep="\t")

#use top1000 and bottom1000 for evaluation
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

data$Training<-rep(FALSE,length(data$Gene))
set.seed(1234)
training_top_index<-sample(which(data$Expressed == "Top1000"),300)
set.seed(1234)
training_bottom_index<-sample(which(data$Expressed == "Bottom1000"),300)
data$Training[training_top_index]<-TRUE
data$Training[training_bottom_index]<-TRUE
training_data<-data.frame(Broad=c(data$Broad.TSS.Coverage.Norm[training_top_index],data$Broad.TSS.Coverage.Norm[training_bottom_index]),Small=c(data$Small.TSS.coverage[training_top_index],data$Small.TSS.coverage[training_bottom_index]))
type<-rep("Unexpressed",length(training_data$Broad))
type[1:300]<-"Expressed"
training_data$Type<-as.factor(type)
svm_model<-svm(Type ~ .,data=training_data,type="C-classification")

data$Prediction<-rep("unknown",length(data$Gene))
test_data_top_index<-which(data$Expressed == "Top1000" & data$Training == FALSE)
test_data_bottom_index<-which(data$Expressed == "Bottom1000" & data$Training == FALSE)
test_data<-data.frame(Broad=c(data$Broad.TSS.Coverage.Norm[test_data_top_index],data$Broad.TSS.Coverage.Norm[test_data_bottom_index]),Small=c(data$Small.TSS.coverage[test_data_top_index],data$Small.TSS.coverage[test_data_bottom_index]))
type<-rep("Unexpressed",length(test_data$Broad))
type[1:length(test_data_top_index)]<-"Expressed"
test_data$Type<-as.factor(type)
prediction<-predict(svm_model,test_data)
predicted_genes<-data.frame(name=c(data$Gene[test_data_top_index],data$Gene[test_data_bottom_index]),prediction=prediction)
data$Prediction[c(test_data_top_index,test_data_bottom_index)]<-as.character(prediction)
write.table(table(pred = prediction, true = test_data$Type),file="8_SVM_Top1000/Prediction.txt",row.names=TRUE,quote=FALSE,sep="\t")

cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
cols[training_top_index]<-"red"
cols[training_bottom_index]<-"green"
png("8_SVM_Top1000/MergedControls_Broad_vs_Small_TSS_coverage_TrainingData.png")
plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
dev.off()

true_expressed<-which(data$Expressed == "Top1000" & data$Training == FALSE & data$Prediction == "Expressed")
true_unexpressed<-which(data$Expressed == "Bottom1000" & data$Training == FALSE & data$Prediction == "Unexpressed")
false_expressed<-which(data$Expressed == "Bottom1000" & data$Training == FALSE & data$Prediction == "Expressed")
false_unexpressed<-which(data$Expressed == "Top1000" & data$Training == FALSE & data$Prediction == "Unexpressed")


cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
cols[false_unexpressed]<-"blue"
cols[false_expressed]<-"blue"
cols[true_expressed]<-"red"
cols[true_unexpressed]<-"green"
png("8_SVM_Top1000/MergedControls_Broad_vs_Small_TSS_coverage_PredictionData.png")
plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
dev.off()


gene_names<-data.frame(data$Gene[true_expressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="8_SVM_Top1000/MergedControls_assigned_in_Top1000_in_expressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

gene_names<-data.frame(data$Gene[true_unexpressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="8_SVM_Top1000/MergedControls_assigned_in_Bottom1000_in_unexpressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

gene_names<-data.frame(data$Gene[false_unexpressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="8_SVM_Top1000/MergedControls_assigned_in_Top1000_in_unexpressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

gene_names<-data.frame(data$Gene[false_expressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="8_SVM_Top1000/MergedControls_assigned_in_Bottom1000_in_expressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")


##############################################################################################################################################################################################
#use Housekeeping and unexpressed as training data
housekeeping<-read.csv("../../../ref/Housekeeping/HK_gene_names.txt",header=FALSE)
unexpressed<-read.csv("../../../ref/FANTOM5/Fantom5_all_lower0.1.txt",header=FALSE)

top_in_expressed_cluster<-0
top_in_unexpressed_cluster<-0
bottom_in_expressed_cluster<-0
bottom_in_unexpressed_cluster<-0
#use data$Expressed vector from above
#data$Expressed<-rep("unknown",length(data$Gene))
for (top in housekeeping$V1) {
    index<-which(data$Gene == top)
    if (length(index) == 0) {
        next
    }
    data$Expressed[index] <- "Housekeeping"
}
for (bottom in unexpressed$V1) {
    index<-which(data$Gene == bottom)
    if (length(index) == 0) {
        next
    }
    data$Expressed[index] <- "Unexpressed"
}

data$Training<-rep(FALSE,length(data$Gene))
set.seed(1234)
training_housekeeping_index<-which(data$Expressed == "Housekeeping")
set.seed(1234)
training_unexpressed_index<-which(data$Expressed == "Unexpressed")
data$Training[training_housekeeping_index]<-TRUE
data$Training[training_unexpressed_index]<-TRUE
training_data<-data.frame(Broad=c(data$Broad.TSS.Coverage.Norm[training_housekeeping_index],data$Broad.TSS.Coverage.Norm[training_unexpressed_index]),Small=c(data$Small.TSS.coverage[training_housekeeping_index],data$Small.TSS.coverage[training_unexpressed_index]))
type<-rep("Unexpressed",length(training_data$Broad))
type[1:length(training_housekeeping_index)]<-"Expressed"
training_data$Type<-as.factor(type)
svm_model<-svm(Type ~ .,data=training_data,type="C-classification")

data$Prediction<-rep("unknown",length(data$Gene))
test_data_top_index<-which(data$Expressed == "Top1000" & data$Training == FALSE)
test_data_bottom_index<-which(data$Expressed == "Bottom1000" & data$Training == FALSE)
test_data<-data.frame(Broad=c(data$Broad.TSS.Coverage.Norm[test_data_top_index],data$Broad.TSS.Coverage.Norm[test_data_bottom_index]),Small=c(data$Small.TSS.coverage[test_data_top_index],data$Small.TSS.coverage[test_data_bottom_index]))
type<-rep("Unexpressed",length(test_data$Broad))
type[1:length(test_data_top_index)]<-"Expressed"
test_data$Type<-as.factor(type)
prediction<-predict(svm_model,test_data)
predicted_genes<-data.frame(name=c(data$Gene[test_data_top_index],data$Gene[test_data_bottom_index]),prediction=prediction)
data$Prediction[c(test_data_top_index,test_data_bottom_index)]<-as.character(prediction)
write.table(table(pred = prediction, true = test_data$Type),file="9_SVM_Housekeeping/Prediction.txt",row.names=TRUE,quote=FALSE,sep="\t")

cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
cols[training_housekeeping_index]<-"red"
cols[training_unexpressed_index]<-"green"
png("9_SVM_Housekeeping/MergedControls_Broad_vs_Small_TSS_coverage_TrainingData.png")
plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
dev.off()

true_expressed<-which(data$Expressed == "Top1000" & data$Training == FALSE & data$Prediction == "Expressed")
true_unexpressed<-which(data$Expressed == "Bottom1000" & data$Training == FALSE & data$Prediction == "Unexpressed")
false_expressed<-which(data$Expressed == "Bottom1000" & data$Training == FALSE & data$Prediction == "Expressed")
false_unexpressed<-which(data$Expressed == "Top1000" & data$Training == FALSE & data$Prediction == "Unexpressed")


cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
cols[false_unexpressed]<-"blue"
cols[false_expressed]<-"blue"
cols[true_expressed]<-"red"
cols[true_unexpressed]<-"green"
png("9_SVM_Housekeeping/MergedControls_Broad_vs_Small_TSS_coverage_PredictionData.png")
plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
dev.off()

cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
cols[false_unexpressed]<-"red"
cols[false_expressed]<-"green"
png("9_SVM_Housekeeping/MergedControls_Broad_vs_Small_TSS_coverage_FalsePrediction.png")
plot(data$Broad.TSS.Coverage.Norm,data$Small.TSS.coverage,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,1.5),col=cols,ylim=c(0,1.5))
dev.off()

gene_names<-data.frame(data$Gene[true_expressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="9_SVM_Housekeeping/MergedControls_assigned_in_Top1000_in_expressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

gene_names<-data.frame(data$Gene[true_unexpressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="9_SVM_Housekeeping/MergedControls_assigned_in_Bottom1000_in_unexpressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

gene_names<-data.frame(data$Gene[false_unexpressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="9_SVM_Housekeeping/MergedControls_assigned_in_Top1000_in_unexpressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")

gene_names<-data.frame(data$Gene[false_expressed], stringsAsFactors=FALSE)
names(gene_names)<-c("GeneName")
write.table(gene_names$GeneName, file="9_SVM_Housekeeping/MergedControls_assigned_in_Bottom1000_in_expressed_cluster.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\n")


