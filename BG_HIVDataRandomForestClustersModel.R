

##Basal Ganglia Analysis File 5
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

setwd("MD/HIVData/ModelAfterRevision")
load("HIVDataTrainTestDataCreation.RData")
load("BasalGnagliaAnalysis.RData")
load ("BG_gxnaTrainClusters.RData")
load ("ModeAfterRevisionGXNAInputCreation.RData")
load("BG_HIVRFClusterModelAfterRevision.RData")



###Do the file creation .phe,.ann. and .exp file creations based on the GXNA tool 
###Save the files at C:\Users\18134\Desktop\gxna
##load the mm vector from the R file called "HIVDataGXNAInputFilesAfterRevision.R" saved at C:\Users\18134\OneDrive - The University of Texas-Rio Grande Valley\Dr.Upal\Bryan\MD\HIVData\ModelAfterRevision

### Load the initial clusters related to Basal Gangliafrom GXNA into an array called "mm_BG[[]]"
load("ModeAfterRevisionGXNAInputCreation.RData")

##Finally we run the RF algorithm. NOTE: we use an ODD number for ntree. This is because when the forest/ensembl is used on test data, ties are broken randomly. 
##Having an odd number of trees avoids this issue and makes the model fully deterministic. 
##Also note, we will use down-sampling to attempt to compensate for unequal class-sizes (less ##relapses than non-relapses).

tmp_BG = as.vector(table(trainY))
num_classes_BG = length(tmp_BG)
min_size_BG = tmp_BG[order(tmp_BG,decreasing=FALSE)[1]]
sampsizes_BG = rep(min_size_BG,num_classes_BG)


install.packages("randomForest")
library(randomForest)

which(filt_genes[,3]=="10379")
which(filt_genes[,3]=="CTGF")
filt_genes[395,1:3]

########################################Model for White Matter

##Part I
##store the Single gene clusters to "singlegene_BG" and multiple gene clusters to "clusters into the"clusterVar_BG". 
singlegene_BG=character(50)
clusterVar_BG<-matrix(list(),1,50)
x=1;
y=1;
for(j in 1:50)
{if (dim(mm_BG[[j]])[1]==1) {
  #singlegene[x]=as.character(mm[[j]]$V2)
  singlegene_BG[x]=mm_BG[[j]]$V2
  x=x+1}
  else
  {clusterVar_BG[[y]]=mm_BG[[j]]$V2
  y=y+1}
}


###Part II
##Creating a matrix to denote the "control", "HIV", "HAD" and "HIVE" in the trainY (we have 19 subjects in the train set)
xx=matrix(0,nrow=19, ncol=4)
xx
for(i in 1:19){j=1
if(trainY[[i]]=="control"){xx[i,j]=1}
else if(trainY[[i]]=="HIV"){xx[i,j+1]=1}
else if(trainY[[i]]=="HAD"){xx[i,j+2]=1}
else if(trainY[[i]]=="HIVE"){xx[i,j+3]=1}
}
####################################
##Part III
###Storing Each gene cluster to "filter_BG" matrix

install.packages("randomForest")
library(randomForest)

filter_BG<-matrix(list(),1,11)
nVar_BG=rep(0,11)
cc_BG<-matrix(list(),11,1)
rf_output_BG<-matrix(list(),11,1)
oberr_BG<-matrix(list(),11,1)
xxx_BG=rep(0,50)
sumxxx_BG=0
#Multiple genes clusters random Forest Models creations
#Clustervar_BG has values upto 11th cell. That is why j=1:11.
##We may have to see how many genes within each cluster and set the following like this"  "else if (nVar_BG[j]==Number of genes)"
for(j in 1:11){
  nVar_BG[j]=length(clusterVar_BG[[j]])
  asw_BG<-as.character(clusterVar_BG[[j]])
  if (nVar_BG[j]==2){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],"Age","Ethnicity"))
  }
  else if (nVar_BG[j]==3){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],asw_BG[3],"Age","Ethnicity"))
  }
  else if (nVar_BG[j]==4){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],asw_BG[3],asw_BG[4],"Age","Ethnicity"))
  }
  else if (nVar_BG[j]==5){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],asw_BG[3],asw_BG[4],asw_BG[5],"Age","Ethnicity"))
  }
  else if (nVar_BG[j]==6){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],asw_BG[3],asw_BG[4],asw_BG[5],asw_BG[6],"Age","Ethnicity"))
  }
  else if (nVar_BG[j]==8){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],asw_BG[3],asw_BG[4],asw_BG[5],asw_BG[6],asw_BG[7],asw_BG[8],"Age","Ethnicity"))
  }
  else if (nVar_BG[j]==15){
    filter_BG[[j]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw_BG[1],asw_BG[2],asw_BG[3],asw_BG[4],asw_BG[5],asw_BG[6],asw_BG[7],asw_BG[8],asw_BG[9],asw_BG[10],asw_BG[11],asw_BG[12],asw_BG[13],asw_BG[14],asw_BG[15],"Age","Ethnicity")) # dim(newpredictor_traindata1GeneID_BG)=19 by 1481
  }
  
  rf_BG=randomForest(x=as.data.frame(filter_BG[[j]]), y=trainY, importance = TRUE, ntree = 10001, proximity=TRUE, na.action = na.omit,mtry=nVar_BG[j])
  for(i in 1:19){
    xxx_BG[[i]]=-(log(rf_BG$votes[[i,1]])*xx[i,1]+log(rf_BG$votes[[i,2]])*xx[i,2]+log(rf_BG$votes[[i,3]])*xx[i,3]+log(rf_BG$votes[[i,4]])*xx[i,4])
    sumxxx_BG=sumxxx_BG+xxx_BG[[i]]
    cc_BG[[j,1]]=sumxxx_BG
    rf_output_BG[[j,1]]=rf_BG
    oberr_BG[[j,1]]=rf_output_BG[[j,1]]$err.rate[10001]
  }
}

##Single Genes Random Forest Model Creations by merging all single genes as one cluster
filter2_BG<-matrix(list(),1,1)
singlegene2_BG=singlegene_BG[1:39]##39 is the number of single gens we had
asw2_BG<-as.character(singlegene2_BG)
filter2_BG[[1]]=subset(newpredictor_traindata1GeneID_BG,select=c(asw2_BG,"Age","Ethnicity"))
nVar2_BG=dim(filter2_BG[[1]])[2]
sumxxx_BG=0

cc2_BG<-matrix(list(),1,nVar2_BG)
rf_output2_BG<-matrix(list(),1,nVar2_BG)
oberr2_BG<-matrix(list(),1,nVar2_BG)

for(mtry in 1:nVar2_BG){
  rf_BG=randomForest(x=as.data.frame(filter2_BG[[1]]), y=trainY, importance = TRUE, ntree = 10001, proximity=TRUE, na.action = na.omit,mtry=mtry)
  for(i in 1:19){
    xxx_BG[[i]]=-(log(rf_BG$votes[[i,1]])*xx[i,1]+log(rf_BG$votes[[i,2]])*xx[i,2]+log(rf_BG$votes[[i,3]])*xx[i,3]+log(rf_BG$votes[[i,4]])*xx[i,4])
    sumxxx_BG=sumxxx_BG+xxx_BG[[i]]
  }
  cc2_BG[[1,mtry]]=sumxxx_BG
  rf_output2_BG[[1,mtry]]=rf_BG
  oberr2_BG[[1,mtry]]=rf_output2_BG[[1,mtry]]$err.rate[10001]
}

####################################
##Part IV
# Transfer error to one error
Finalerror_BG=rep(0,12) ##23 clusters + single gene cluster
for(i in 1:12){
  if (i<12){
    Finalerror_BG[i]=oberr_BG[[i,1]]}
  else 
  { Finalerror_BG[i]=oberr2_BG[[1,5]]}##Oberr2_BG[[1,5]] has the lowest error. That is why mtry=5 has been selected.
}
Finalerror_BG

######################## transfer data clusters to one data frame
DataFinal_BG<-matrix(list(),1,12)
ID_BG=as.vector(rownames(newpredictor_traindata1GeneID_BG))
for(i in 1:24){if(i<12){
  DataFinal_BG[[i]]=filter_BG[[i]]
  DataFinal_BG[[i]]<-cbind(DataFinal_BG[[i]],ID_BG)}
  else {DataFinal_BG[[i]]=filter2_BG[[1]]
  DataFinal_BG[[i]]<-cbind(DataFinal_BG[[i]],ID_BG)}
}

#### Look at Cluster OOB Error 
err.temp_BG <- abs(Finalerror_BG) 
err.temp_BG
errt_BG <- order(err.temp_BG,decreasing=FALSE) 
errt_BG
colnames(DataFinal_BG[[12]])  ##Cluster Ranked Number 1. However, this is the single gene cluster model and Hence we do not use it.
colnames(DataFinal_BG[[3]])  ##Cluster Ranked Number 2
colnames(DataFinal_BG[[4]])  ##Cluster Ranked Number 3
colnames(DataFinal_BG[[10]])  ##Cluster Ranked Number 4
colnames(DataFinal_BG[[6]])  ##Cluster Ranked Number 5


dev.off()
plot(c(1:length(Finalerror_BG)),err.temp_BG[errt_BG],log='x',cex.main=1.5, xlab='Cluster Rank',ylab='Cluster Importance (OOB Error)',cex.lab=1, pch=16,main='Cluster Importance for HIV Using RF') 

errt_BG
rf_output_BG[[j,1]]

rf_importances1_BG=importance(rf_output2_BG[[1,5]], scale=FALSE) ##This is for the Singke gene Cluster Model. That is there is only 1 row. coulmn is the best mtry value with the lowest error
rf_importances2_BG=importance(rf_output_BG[[3]], scale=FALSE)
rf_importances3_BG=importance(rf_output_BG[[4]], scale=FALSE)
rf_importances4_BG=importance(rf_output_BG[[10]], scale=FALSE)

# Predicting response variable
predicted.response1_BG <- predict(rf_output2_BG[[1,5]],DataFinal_BG[[12]])
table(trainY,predicted.response1_BG)
predicted.response2_BG <- predict(rf_output_BG[[3]],DataFinal_BG[[3]])
table(trainY,predicted.response2_BG)
predicted.response3_BG <- predict(rf_output_BG[[4]],DataFinal_BG[[4]])
table(trainY,predicted.response3_BG)
predicted.response4_BG <- predict(rf_output_BG[[10]],DataFinal_BG[[10]])
table(trainY,predicted.response4_BG)

###This has been used in this study
predicted.testresponse1_BG <- predict(rf_output2_BG[[1,5]],newpredictor_testdata1GeneID_BG)
table(testY,predicted.testresponse1_BG)
aa1_BG=table(predicted.testresponse1_BG,testY)
Accuracy1_BG=(aa1_BG[1,1]+aa1_BG[2,2]+aa1_BG[3,3]+aa1_BG[4,4])/sum(aa1_BG)
Accuracy1_BG
predicted.testresponse2_BG <- predict(rf_output_BG[[3]],newpredictor_testdata1GeneID_BG)
table(testY,predicted.testresponse2_BG)
aa2_BG=table(predicted.testresponse2_BG,testY)
Accuracy2_BG=(aa2_BG[1,1]+aa2_BG[2,2]+aa2_BG[3,3]+aa2_BG[4,4])/sum(aa2_BG)
Accuracy2_BG
predicted.testresponse3_BG <- predict(rf_output_BG[[4]],newpredictor_testdata1GeneID_BG)
table(testY,predicted.testresponse3_BG)
aa3_BG=table(predicted.testresponse3_BG,testY)
Accuracy3_BG=(aa3_BG[1,1]+aa3_BG[2,2]+aa3_BG[3,3]+aa3_BG[4,4])/sum(aa3_BG)
Accuracy3_BG##So far, the best Accuracy
predicted.testresponse4_BG <- predict(rf_output_BG[[10]],newpredictor_testdata1GeneID_BG)
table(testY,predicted.testresponse4_BG)
aa4_BG=table(predicted.testresponse4_BG,testY)
Accuracy4_BG=(aa4_BG[1,1]+aa4_BG[2,2]+aa4_BG[3,3]+aa4_BG[4,4])/sum(aa4_BG)
Accuracy4_BG

plot.new()
plot(1:SinglenVar1,singleoberr1_BG,pch=19,col=c("red"),type="b",ylim=c(0,4),ylab="OOB Error",xlab="Number of Predictors Considered at each Split")
legend("topright",legend=c("Out of Bag Error"),pch=19, col=c("red"))

###Using the best Accuracy model
pred_BG<-predict(rf_output_BG[[4]],newpredictor_testdata1GeneID_BG)
aa_BG=table(pred_BG,testY)
aa_BG
pred_BG2<-predict(rf_output_BG[[3]],newpredictor_testdata1GeneID_BG)
aa_BG2=table(pred_BG,testY)
aa_BG2

###Model Evaluation Measures
##Accuracy is another measure that can be useful when the problem has well balanced classes (for example in optical character recognition) and we want to put an emphasis on exact matches. 
##Unforunately accuracy suffers on unbalanced data sets, a typical example is information retrieval.
Accuracy_BG=(aa_BG[1,1]+aa_BG[2,2]+aa_BG[3,3]+aa_BG[4,4])/sum(aa_BG)
Accuracy_BG
Accuracy_BG2=(aa_BG2[1,1]+aa_BG2[2,2]+aa_BG2[3,3]+aa_BG2[4,4])/sum(aa_BG2)
Accuracy_BG2

###Precision:Intuitively, a high precision for a class means that if our models predict that class, it is very likely to be true. 
###A high precision model will be useful in those situations where we need to have an high confidence in our prediction (for example in medical diagnosis).
pre_c_BG=aa_BG[1,1]/sum(aa_BG[1,])
pre_HAD_BG=aa_BG[2,2]/sum(aa_BG[2,])
pre_HIV_BG=aa_BG[3,3]/sum(aa_BG[3,])
pre_HIVE_BG=aa_BG[4,4]/sum(aa_BG[4,])
pre_c_BG
pre_HAD_BG
pre_HIV_BG
pre_HIVE_BG

pre_c_BG2=aa_BG2[1,1]/sum(aa_BG2[1,])
pre_HAD_BG2=aa_BG2[2,2]/sum(aa_BG2[2,])
pre_HIV_BG2=aa_BG2[3,3]/sum(aa_BG2[3,])
pre_HIVE_BG2=aa_BG2[4,4]/sum(aa_BG2[4,])
pre_c_BG2
pre_HAD_BG2
pre_HIV_BG2
pre_HIVE_BG2
###REcall: If recall is high, it means that our models manages to recover most instances of that class. 
Re_c_BG=aa_BG[1,1]/sum(aa_BG[,1])
Re_HAD_BG=aa_BG[2,2]/sum(aa_BG[,2])
Re_HIV_BG=aa_BG[3,3]/sum(aa_BG[,3])
Re_HIVE_BG=aa_BG[4,4]/sum(aa_BG[,4])
Re_c_BG
Re_HAD_BG
Re_HIV_BG
Re_HIVE_BG

Re_c_BG2=aa_BG2[1,1]/sum(aa_BG2[,1])
Re_HAD_BG2=aa_BG2[2,2]/sum(aa_BG2[,2])
Re_HIV_BG2=aa_BG2[3,3]/sum(aa_BG2[,3])
Re_HIVE_BG2=aa_BG2[4,4]/sum(aa_BG2[,4])
Re_c_BG2
Re_HAD_BG2
Re_HIV_BG2
Re_HIVE_BG2
###F1-score
##F1 score is the harmonic mean of precision and recall, and acts as a combined measure of the two.

F1_c_BG=(2*pre_c_BG*Re_c_BG)/(pre_c_BG+Re_c_BG)
F1_HAD_BG=(2*pre_HAD_BG*Re_HAD_BG)/(pre_HAD_BG+Re_HAD_BG)
F1_HIV_BG=(2*pre_HIV_BG*Re_HIV_BG)/(pre_HIV_BG+Re_HIV_BG)
F1_HIVE_BG=(2*pre_HIVE_BG*Re_HIVE_BG)/(pre_HIVE_BG+Re_HIVE_BG)


###Macro Average
macPrecision_BG=(pre_c_BG+pre_HAD_BG+pre_HIV_BG+pre_HIVE_BG)/4
macPrecision_BG
macPrecision_BG=(pre_c_BG+pre_HIV_BG+pre_HIVE_BG)/3 ##Since pre_HAD_BG=NaN
macPrecision_BG
macRecall_BG=(Re_c_BG+Re_HAD_BG+Re_HIV_BG+Re_HIVE_BG)/4
macRecall_BG

macPrecision_BG2=(pre_c_BG2+pre_HAD_BG2+pre_HIV_BG2+pre_HIVE_BG2)/4
macPrecision_BG2
macPrecision_BG2=(pre_HAD_BG2+pre_HIVE_BG2)/2 ##Since pre_HIV_FC=NaN and pre_c_FC=NaN
macPrecision_BG2
macRecall_BG2=(Re_c_BG2+Re_HAD_BG2+Re_HIV_BG2+Re_HIVE_BG2)/4
macRecall_BG2

###Micro Average
micPrecision_BG=(pre_c_BG+pre_HAD_BG+pre_HIV_BG+pre_HIVE_BG)/sum(aa_BG[,1]+aa_BG[,2]+aa_BG[,3]+aa_BG[,4])
micPrecision_BG

micPrecision_BG=(pre_c_BG+pre_HAD_BG+pre_HIVE_BG)/sum(aa_BG[,1]+aa_BG[,2]+aa_BG[,3]+aa_BG[,4]) ##Since pre_HIV_BG=NaN
micPrecision_BG

micRecall_BG=(Re_c_BG+Re_HAD_BG+Re_HIV_BG+Re_HIVE_BG)/sum(aa_BG[,1]+aa_BG[,2]+aa_BG[,3]+aa_BG[,4])
micRecall_BG


# Look at Gene (Variables) importance 


##This has been used in this program
##Using Top 2 and3 clusters
imp.tempsingle3_BG <- abs(rf_output_BG[[3]]$importance[,5]) 
write.table(imp.tempsingle3_BG,"imp.tempsingle3_BG.txt")
single3t_BG <- order(imp.tempsingle3_BG,decreasing=TRUE) 

imp.tempsingle4_BG <- abs(rf_output_BG[[4]]$importance[,5]) 
write.table(imp.tempsingle4_BG,"imp.tempsingle4_BG.txt")
single4t_BG <- order(imp.tempsingle4_BG,decreasing=TRUE) 


##Double Checking Gene ID vs Gene Name
geneAnnot_BG[which(geneAnnot_BG[,2]==3959),1]

dev.off()
##DataFinal_BG[[3]][,-11])) remove column 11 as it is "ID_BG" variable. This has 8 gens +Age+Ethnicity=Toal 10 variables. 
plot(c(1:ncol(DataFinal_BG[[3]][,-11])),imp.tempsingle3_BG[single3t_BG],log='x',cex.main=1.5, xlab='Gene Rank Using RF Clusters',ylab='Variable Importance (Mean Decrease Accuracy)',cex.lab=1, pch=16,main='Basal Ganglia') 

# Get subset of expression values for 15 most 'important' genes 
gn.imp3_BG <- names(imp.tempsingle3_BG)[single3t_BG] 
Bt_BG<-which(gn.imp3_BG=="Ethnicity")
At_BG<-which(gn.imp3_BG=="Age")
#cumDataNew1<-cumData1[,!(names(cumData1) %in% c("ID"))]
gn.8_BG <- gn.imp3_BG[-c(Bt_BG,At_BG)][1:8] # vector of top 8 genes (we only have 8 genes in this cluster), without Ethnicity and Age
t3_BG <- is.element(colnames(DataFinal_BG[[3]]),gn.8_BG) 
#t1 <- is.element(colnames(cumDataNew1),gn.15) 
#sig.eset1 <- cumDataNew1[,t1]
sig.eset3_BG<- DataFinal_BG[[3]][,t3_BG]
sig.esetNew3_BG<-t(sig.eset3_BG) 

nameArr3_BG<-character(8)
for(j in 1:8){
  nameArr3_BG[j] =geneAnnot_BG[geneAnnot_BG[,2]==rownames(sig.esetNew3_BG)[j],1]
}
nameArr3_BG

#########################
# Get subset of expression values for 15 most 'important' genes 
gn.imp4_BG <- names(imp.tempsingle4_BG)[single4t_BG] 
Bt_BG<-which(gn.imp4_BG=="Ethnicity")
At_BG<-which(gn.imp4_BG=="Age")
#cumDataNew1<-cumData1[,!(names(cumData1) %in% c("ID"))]
gn4.4_BG <- gn.imp4_BG[-c(Bt_BG,At_BG)][1:4] # vector of top 4 genes, without BType and Age
t1_4_BG <- is.element(colnames(DataFinal_BG[[4]]),gn4.4_BG) 
#t1 <- is.element(colnames(cumDataNew1),gn.15) 
#sig.eset1 <- cumDataNew1[,t1]
sig.eset4_BG <- DataFinal_BG[[4]][,t1_4_BG]
sig.esetNew4_BG<-t(sig.eset4_BG) 

nameArr4_BG<-character(4)
for(j in 1:4){
  nameArr4_BG[j] =geneAnnot_BG[geneAnnot_BG[,2]==rownames(sig.esetNew4_BG)[j],1]
}
nameArr4_BG

##########################################Heat Map Creation for All Data##
# matrix of expression values, not necessarily in order ## Make a heatmap, with group differences obvious on plot 
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256) 
rownames(filteredDataAll) <- target # This will label the heatmap columns 

csc <- rep(hmcol[50],24) 
csc[rownames(filteredDataAll)=='HAD'] <- hmcol[200] # column side colors
csc[rownames(filteredDataAll)=='HIV'] <- hmcol[210]
csc[rownames(filteredDataAll)=='HIVE'] <- hmcol[220]

colnames(filteredDataAll_BG)=geneAnnot_BG[,1] ##Naming filteredDataAll columns by gene names

install.packages("gplots")
library(gplots)

filteredDataAllTopGenes3_BG=filteredDataAll_BG[,nameArr3_BG]
dim(filteredDataAllTopGenes3_BG)
filteredDataAllTopGenes4_BG=filteredDataAll_BG[,nameArr4_BG]
dim(filteredDataAllTopGenes4_BG)

dev.off()
heatmap.2(t(filteredDataAllTopGenes3_BG),scale="row",col=hmcol,ColSideColors=csc[1:24],main="Top 15 Gene Expression Variations for Basal Ganglia using RF Cluster 1",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')
dev.off()
heatmap.2(t(filteredDataAllTopGenes4_BG),scale="row",col=hmcol,ColSideColors=csc[1:24],main="Top 15 Gene Expression Variations for Basal Ganglia using RF Cluster 2",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


#########################################

save.image("BG_HIVRFClusterModelAfterRevision.RData")
################