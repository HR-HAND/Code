
##This is the Basal Ganglia File 3

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())


#We need to load the following data from File 1. 
load("HIVDataTrainTestDataCreation.RData")
load("BasalGangliaAnalysis.RData")
load("BG_HIVSingleGenesModelAfterRevision.RData")
install.packages("randomForest")
library(randomForest)

tmp = as.vector(table(trainY))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)
#rf_output$mtry=41 
set.seed(620)
rf_outputSingleGenes_BG=randomForest(x=newpredictor_traindata1_BG, y=trainY, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes, na.action = na.omit)
rf_outputSingleGenes_BG$err.rate[10001,]

dev.off()
plot(rf_outputSingleGenes_BG) ##Black line: OOB error rate ##Other lines: Refres to the error rate for each class
rf_outputSingleGenes_BG$confusion
rf_outputSingleGenes_BG$votes

save(rf_outputSingleGenes_BG, file="RFSingleGenes_BG_model")
load("RFSingleGenes_BG_model")
rf_importancesSingleGene_BG=importance(rf_outputSingleGenes_BG, scale=FALSE)

# Predicting response variable
predicted.responseSingleGenes_BG <- predict(rf_outputSingleGenes_BG,newpredictor_traindata1_BG)
table(trainY,predicted.responseSingleGenes_BG)
predicted.testresponseSingleGenes_BG <- predict(rf_outputSingleGenes_BG,newpredictor_testdata1_BG)
table(testY,predicted.testresponseSingleGenes_BG)


#aa=table(pred, testY)
SingleGenesaa_BG=table(predicted.testresponseSingleGenes_BG,testY)

###Model Evaluation Measures

##Accuracy is another measure that can be useful when the problem has well balanced classes (for example in optical character recognition) and we want to put an emphasis on exact matches. 
##Unforunately accuracy suffers on unbalanced data sets, a typical example is information retrieval.
SingleGenesAccuracy_BG=(SingleGenesaa_BG[1,1]+SingleGenesaa_BG[2,2]+SingleGenesaa_BG[3,3]+SingleGenesaa_BG[4,4])/sum(SingleGenesaa_BG)
SingleGenesAccuracy_BG

###Precision:Intuitively, a high precision for a class means that if our models predict that class, it is very likely to be true. 
###A high precision model will be useful in those situations where we need to have an high confidence in our prediction (for example in medical diagnosis).
SingleGenespre_c_BG=SingleGenesaa_BG[1,1]/sum(SingleGenesaa_BG[1,])
SingleGenespre_HAD_BG=SingleGenesaa_BG[2,2]/sum(SingleGenesaa_BG[2,])
SingleGenespre_HIV_BG=SingleGenesaa_BG[3,3]/sum(SingleGenesaa_BG[3,])
SingleGenespre_HIVE_BG=SingleGenesaa_BG[4,4]/sum(SingleGenesaa_BG[4,])
SingleGenespre_c_BG
SingleGenespre_HAD_BG
SingleGenespre_HIV_BG
SingleGenespre_HIVE_BG

###REcall: If recall is high, it means that our models manages to recover most instances of that class. 
SingleGenesRe_c_BG=SingleGenesaa_BG[1,1]/sum(SingleGenesaa_BG[,1])
SingleGenesRe_HAD_BG=SingleGenesaa_BG[2,2]/sum(SingleGenesaa_BG[,2])
SingleGenesRe_HIV_BG=SingleGenesaa_BG[3,3]/sum(SingleGenesaa_BG[,3])
SingleGenesRe_HIVE_BG=SingleGenesaa_BG[4,4]/sum(SingleGenesaa_BG[,4])
SingleGenesRe_c_BG
SingleGenesRe_HAD_BG
SingleGenesRe_HIV_BG
SingleGenesRe_HIVE_BG

###F1-score
##F1 score is the harmonic mean of precision and recall, and acts as a combined measure of the two.

SingleGenesF1_c_BG=(2*SingleGenespre_c_BG*SingleGenesRe_c_BG)/(SingleGenespre_c_BG+SingleGenesRe_c_BG)
SingleGenesF1_HAD_BG=(2*SingleGenespre_HAD_BG*SingleGenesRe_HAD_BG)/(SingleGenespre_HAD_BG+SingleGenesRe_HAD_BG)
SingleGenesF1_HIV_BG=(2*SingleGenespre_HIV_BG*SingleGenesRe_HIV_BG)/(SingleGenespre_HIV_BG+SingleGenesRe_HIV_BG)
SingleGenesF1_HIVE_BG=(2*SingleGenespre_HIVE_BG*SingleGenesRe_HIVE_BG)/(SingleGenespre_HIVE_BG+SingleGenesRe_HIVE_BG)
SingleGenesF1_c_BG
SingleGenesF1_HAD_BG
SingleGenesF1_HIV_BG
SingleGenesF1_HIVE_BG
###Macro Average
SingleGenesmacPrecision_BG=(SingleGenespre_c_BG+SingleGenespre_HAD_BG+SingleGenespre_HIV_BG+SingleGenespre_HIVE_BG)/4  ##Since +SingleGenespre_HIV_BG=NaN
SingleGenesmacRecall_BG=(SingleGenesRe_c_BG+SingleGenesRe_HAD_BG+SingleGenesRe_HIV_BG+SingleGenesRe_HIVE_BG)/4
SingleGenesmacPrecision_BG
SingleGenesmacRecall_BG

###Micro Average
SingleGenesmicPrecision_BG=(SingleGenespre_c_BG+SingleGenespre_HAD_BG+SingleGenespre_HIV_BG+SingleGenespre_HIVE_BG)/sum(SingleGenesaa_BG[,1]+SingleGenesaa_BG[,2]+SingleGenesaa_BG[,3]+SingleGenesaa_BG[,4])
SingleGenesmicRecall_BG=(SingleGenesRe_c_BG+SingleGenesRe_HAD_BG+SingleGenesRe_HIV_BG+SingleGenesRe_HIVE_BG)/sum(SingleGenesaa_BG[,1]+SingleGenesaa_BG[,2]+SingleGenesaa_BG[,3]+SingleGenesaa_BG[,4])
SingleGenesmacPrecision_BG
SingleGenesmacRecall_BG
SingleGenesmicPrecision_BG
SingleGenesmicRecall_BG


# Look at variable importance 
SingleGenes_imp.temp_BG <- abs(rf_outputSingleGenes_BG$importance[,5]) 
SingleGenes_t_BG <- order(SingleGenes_imp.temp_BG,decreasing=TRUE) 

dev.off()
plot(c(1:ncol(newpredictor_traindata1_BG)),SingleGenes_imp.temp_BG[SingleGenes_t_BG],log='x',cex.main=1.5, xlab='Gene Rank',ylab='Variable Importance (Mean Decrease Accuracy)',cex.lab=1, pch=16,main='Basal Ganglia Single Genes RF Model') 

# Get subset of expression values for 25 most 'important' genes 
SingleGenes_gn.imp_BG <- names(SingleGenes_imp.temp_BG)[SingleGenes_t_BG] 
SingleGenes_gn.25_BG <- SingleGenes_gn.imp_BG[1:25] # vector of top 25 genes, in order 
SingleGenes_t_BG <- is.element(colnames(newpredictor_traindata1_BG),SingleGenes_gn.25_BG)

SingleGenes_sig.eset_BG <- newpredictor_traindata1_BG[,SingleGenes_t_BG]
SingleGenes_sig.esetNew_BG<-t(SingleGenes_sig.eset_BG)##Dim=25 by 19
str(SingleGenes_sig.esetNew_BG)
dim(SingleGenes_sig.esetNew_BG)

# matrix of expression values, not necessarily in order ## Make a heatmap, with group differences obvious on plot 
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256) 
colnames(SingleGenes_sig.esetNew_BG) <- trainY # This will label the heatmap columns 
rownames(SingleGenes_sig.esetNew_BG)
csc <- rep(hmcol[50],19) 
csc[colnames(SingleGenes_sig.esetNew_BG)=='HAD'] <- hmcol[200] # column side color will be purple for T and orange for B 
csc[colnames(SingleGenes_sig.esetNew_BG)=='HIV'] <- hmcol[210]
csc[colnames(SingleGenes_sig.esetNew_BG)=='HIVE'] <- hmcol[220]


install.packages("gplots")
library(gplots)


dev.off()
heatmap.2(SingleGenes_sig.esetNew_BG,scale="row",col=hmcol,ColSideColors=csc,main="Top 25 Genes in Basal Ganglia Using Single Genes RF",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


##########################################Heat Map Creation for All Data##
# matrix of expression values, not necessarily in order ## Make a heatmap, with group differences obvious on plot 

colnames(filteredDataAll_BG)=predictor_names_BG

library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256) 
rownames(filteredDataAll_BG) <- target # This will label the heatmap columns 

csc <- rep(hmcol[50],24) 
csc[rownames(filteredDataAll_BG)=='HAD'] <- hmcol[200] # column side color will be purple for T and orange for B 
csc[rownames(filteredDataAll_BG)=='HIV'] <- hmcol[210]
csc[rownames(filteredDataAll_BG)=='HIVE'] <- hmcol[220]


install.packages("gplots")
library(gplots)

filteredDataAllTopGenesRFSingle_BG=filteredDataAll_BG[,rownames(SingleGenes_sig.esetNew_BG)]
dim(filteredDataAllTopGenesRFSingle_BG)


dev.off()
heatmap.2(t(filteredDataAllTopGenesRFSingle_BG),scale="row",col=hmcol,ColSideColors=csc[1:24],main="Top 25 Gene Expression Variations for Basal Ganglia",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


save.image("BG_HIVSingleGenesModelAfterRevision.RData")

########################################################
##

# Look at variable importance 
imp.temp <- abs(rf_output$importance[,5]) 
t <- order(imp.temp,decreasing=TRUE) 

# Get subset of expression values for 25 most 'important' genes 
gn.imp <- names(imp.temp)[t] 
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order 
t <- is.element(colnames(newpredictor_traindata1),gn.25)

sig.eset <- newpredictor_traindata1[,t]
sig.esetNew<-t(sig.eset)##Dim=25 by 52
str(sig.esetNew)



# Get subset of expression values for 25 most 'important' genes 
gn.imp1 <- names(imp.tempsingle1)[single1t] 
Bt<-which(gn.imp1=="BType")
At<-which(gn.imp1=="Age")
#cumDataNew1<-cumData1[,!(names(cumData1) %in% c("ID"))]
gn.15 <- gn.imp1[-c(Bt,At)][1:15] # vector of top 15 genes, without BType and Age
t1 <- is.element(colnames(DataFinal[[7]]),gn.15) 
#t1 <- is.element(colnames(cumDataNew1),gn.15) 
#sig.eset1 <- cumDataNew1[,t1]
sig.eset1 <- DataFinal[[7]][,t1]
sig.esetNew1<-t(sig.eset1) 

nameArr1<-character(15)
for(j in 1:15){
  nameArr1[j] =geneAnnot[geneAnnot[,2]==rownames(sig.esetNew1)[j],1]
}

#########################
# Get subset of expression values for 25 most 'important' genes 
gn.imp3 <- names(imp.tempsingle3)[single3t] 
Bt<-which(gn.imp3=="BType")
At<-which(gn.imp3=="Age")
#cumDataNew1<-cumData1[,!(names(cumData1) %in% c("ID"))]
gn3.15 <- gn.imp3[-c(Bt,At)][1:15] # vector of top 15 genes, without BType and Age
t1_3 <- is.element(colnames(DataFinal[[6]]),gn3.15) 
#t1 <- is.element(colnames(cumDataNew1),gn.15) 
#sig.eset1 <- cumDataNew1[,t1]
sig.eset3 <- DataFinal[[6]][,t1_3]
sig.esetNew3<-t(sig.eset3) 

nameArr3<-character(15)
for(j in 1:15){
  nameArr3[j] =geneAnnot[geneAnnot[,2]==rownames(sig.esetNew3)[j],1]
}
##########################################Heat Map Creation for All Data##
# matrix of expression values, not necessarily in order ## Make a heatmap, with group differences obvious on plot 
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256) 
rownames(filteredDataAll) <- target # This will label the heatmap columns 

csc <- rep(hmcol[50],52) 
csc[colnames(filteredDataAll)=='HAD'] <- hmcol[200] # column side color will be purple for T and orange for B 
csc[colnames(filteredDataAll)=='HIV'] <- hmcol[210]
csc[colnames(filteredDataAll)=='HIVE'] <- hmcol[220]


install.packages("gplots")
library(gplots)

filteredDataAllTopGenes2=filteredDataAll[,nameArr2]
dim(filteredDataAllTopGenes2)
filteredDataAllTopGenes3=filteredDataAll[,nameArr3]
dim(filteredDataAllTopGenes3)

dev.off()
heatmap.2(t(filteredDataAllTopGenes2[1:24,]),scale="row",col=hmcol,ColSideColors=csc[1:24],main="Top 15 Gene Expression Variations for White Matter using 1 RF Cluster",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')
dev.off()
heatmap.2(t(filteredDataAllTopGenes3[1:24,]),scale="row",col=hmcol,ColSideColors=csc[1:24],main="Top 15 Gene Expression Variations for White Matter using 1 RF Cluster",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


dev.off()
heatmap.2(t(filteredDataAllTopGenes2[25:48,]),scale="row",col=hmcol,ColSideColors=csc[25:48],main="Top 15 Gene Expression Variations for Frontal Cortex using 1 RF Cluster",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')
dev.off()
heatmap.2(t(filteredDataAllTopGenes3[25:48,]),scale="row",col=hmcol,ColSideColors=csc[25:48],main="Top 15 Gene Expression Variations for Frontal Cortex using 1 RF Cluster",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


dev.off()
heatmap.2(t(filteredDataAllTopGenes2[49:72,]),scale="row",col=hmcol,ColSideColors=csc[49:72],main="Top 15 Gene Expression Variations for Basal Ganglia using 1 RF Cluster",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')
dev.off()
heatmap.2(t(filteredDataAllTopGenes3[49:72,]),scale="row",col=hmcol,ColSideColors=csc[49:72],main="Top 15 Gene Expression Variations for Basal Ganglia using 1 RF Cluster",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


save.image("HIVDataNewRFModel.RData")

