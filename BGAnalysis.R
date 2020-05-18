##DataFile: file:///C:/Users/18134/OneDrive%20-%20The%20University%20of%20Texas-Rio%20Grande%20Valley/data2/GEO%20DataSet%20Browser.html

##ftp://ftp.ncbi.nlm.nih.gov/geo/

###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35864 ##HIV Data
##https://nntc.org/content/gene_array
###SAM Analysis: https://mdozmorov.github.io/BIOS567/assets/presentation_diffexpression/DiffExpr_SAM.html
##Gene Filtering http://math.usu.edu/jrstevens/stat5570/3.2.Filtering.pdf

###Basal Gnaglia Analysis
##File 2

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

setwd("MD/HIVData/ModelAfterRevision")
load("HIVDataTrainTestDataCreation.RData")

load("BasalGangliaAnalysis.RData")



#Get potential predictor variables for Basal Ganglia
newpredictor_traindata_BG=t(filt_TrainData_BG)  ##dim(newpredictor_traindata_BG)=19 by 710
newpredictor_testdata_BG=t(filt_TestData_BG)  ##dim(newpredictor_testdata_BG)=5 by 710

#TrainXDem[1:3,]
#age.ch1 race.ch1   tissue.ch1
#2      53    White White matter
#3      63 Hispanic White matter
#4      58    White White matter

#TrainXDem[,1]=age.ch1
#TrainXDem[,2]=race.ch1
#is.numeric(TrainXDem$age.ch1)
#TRUE

newpredictor_traindata1_BG=data.frame(newpredictor_traindata_BG,as.numeric(TrainXDem[,1]),TrainXDem[,2],stringsAsFactors=FALSE)
newpredictor_testdata1_BG=data.frame(newpredictor_testdata_BG,as.numeric(TestXDem[,1]),TestXDem[,2],stringsAsFactors=FALSE)
##dim(newpredictor_traindata1_BG)= 19 by 712
##dim(newpredictor_testdata1_BG)= 5 by 712

newpredictor_traindata1GeneID_BG=data.frame(newpredictor_traindata_BG,as.numeric(TrainXDem[,1]),TrainXDem[,2], stringsAsFactors=FALSE)
newpredictor_testdata1GeneID_BG=data.frame(newpredictor_testdata_BG,as.numeric(TestXDem[,1]),TestXDem[,2],stringsAsFactors=FALSE)

############################
predictor_dataCombined_BG=data.frame(predictor_dataBG,clincData24Individual[,3],target, stringsAsFactors=FALSE)
colnames(predictor_dataCombined_BG)=c(colnames(predictor_dataCombined_BG[,-c(20409:20410)]),"Age","Target")
##20409:20410 are the varoavbels "Age" and "Target". In order to rename 20409 by Age and 20410 by Target we do this.
## colnames(predictor_dataCombined_BG[1:6,20405:20410])
##"X20416" "X20417" "X20418" "X20419" "Age"    "Target"
##########################################

#Get potential predictor variables


predictor_names_BG=c(as.vector(filt_genes_BG[,3])) #gene symbol for BG
predictor_namesID_BG=c(as.vector(filt_genes_BG[,2])) #gene ID for BG
geneAnnot_BG=cbind(predictor_names_BG,predictor_namesID_BG)

colnames(newpredictor_traindata1GeneID_BG)=c(predictor_namesID_BG,"Age","Ethnicity")
colnames(newpredictor_testdata1GeneID_BG)=c(predictor_namesID_BG,"Age","Ethnicity")

colnames(newpredictor_traindata1_BG)=c(predictor_names_BG,"Age","Ethnicity")
colnames(newpredictor_testdata1_BG)=c(predictor_names_BG,"Age","Ethnicity")


#> newpredictor_traindata1GeneID_BG[1:5,708:712]
#                         9924     9935      998 Age Ethnicity
#GSM876887_50G759.CEL 5.425305 4.232665 5.380421  53     White
#GSM876888_51G759.CEL 6.044727 3.208320 6.334287  63  Hispanic
#GSM876889_52G759.CEL 6.542343 4.267008 5.262936  58     White
#GSM876890_53G759.CEL 5.054077 2.697737 6.394291  34     White
#GSM876891_54G759.CEL 5.054077 4.606351 6.721265  48     White

#> newpredictor_traindata1_BG[1:5,708:712]
#                         PAN2     MAFB    CDC42 Age Ethnicity
#GSM876887_50G759.CEL 5.425305 4.232665 5.380421  53     White
#GSM876888_51G759.CEL 6.044727 3.208320 6.334287  63  Hispanic
#GSM876889_52G759.CEL 6.542343 4.267008 5.262936  58     White
#GSM876890_53G759.CEL 5.054077 2.697737 6.394291  34     White
#GSM876891_54G759.CEL 5.054077 4.606351 6.721265  48     White


dim(newpredictor_traindata1_BG) #19 by 712
dim(newpredictor_testdata1_BG)  #5 by 712
str(newpredictor_testdata1_BG[,708:712])##Checking whether Age variable is taken as numeric and Ethnicity as a fcator

##Testing whether the column names of testing data set is same as for training data set
which(colnames(newpredictor_testdata1_BG)!=colnames(newpredictor_traindata1_BG))
str(colnames(newpredictor_testdata1_BG[712]))

###Creating a combined dataset for BG with both training and Testing data sets
newpredictor_trainTestcombined_BG=rbind(newpredictor_traindata1_BG,newpredictor_testdata1_BG)
newpredictor_trainTestcombined_BG_SAMAnalysis=cbind(as.data.frame.array(t(colnames(newpredictor_traindata1GeneID_BG)),t(newpredictor_trainTestcombined_BG)))


save(newpredictor_traindata1_BG, file="DataFilesnewpredictor_traindata1_BG")
write.table(newpredictor_traindata1_BG, "newpredictor_traindata1_BG.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_traindata1GeneID_BG, file="newpredictor_traindata1GeneID_BG")
write.table(newpredictor_traindata1GeneID_BG, "newpredictor_traindata1GeneID_BG.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_testdata1_BG, file="newpredictor_testdata1_BG")
write.table(newpredictor_testdata1_BG, "newpredictor_testdata1_BG.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_testdata1GeneID_BG, file="newpredictor_testdata1GeneID_BG")
write.table(newpredictor_testdata1GeneID_BG, "newpredictor_testdata1GeneID_BG.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_trainTestcombined_BG, file="newpredictor_trainTestcombined_BG")
write.table(newpredictor_trainTestcombined_BG,"newpredictor_trainTestcombined_BG.txt", quote = FALSE, sep = "\t", row.names = FALSE)


save.image("HIVDataTrainTestDataCreation.RData")

###################################################################################################
setwd("ModelAfterRevision")
load("HIVDataTrainTestDataCreation.RData")

###Multidimensional scaling plot of distances between gene expression profiles
BiocManager::install(version = "3.9")
BiocManager::install("limma")
library(limma)

namessamples_BG<-seq(838,861,by=1) ##Naming by GEO Accession Number
rownames(filteredDataAll_BG)=namessamples_BG

plotMDS(t(filteredDataAll_BG),col=c(rep("black",6), rep("orange",6),rep("purple",7),rep("red",5)), labels= c(rep("Control",6), rep("HIV",6),rep("HAD",7),rep("HIVE",5)))

###Log 2fold Calculation

targetAll<-as.factor(c(trainY,testY))
newpredictor_trainTestcombinedAll_BG=cbind(newpredictor_trainTestcombined_BG,targetAll)
colnames(newpredictor_trainTestcombinedAll_BG)=c(colnames(newpredictor_traindata1_BG),"target")

geneID_BG=c(predictor_namesID_BG)
genenames_BG=colnames(newpredictor_traindata1_BG)[-c(711,712)] ##Filtering only the gene names withoot columns 1482:Age, 1483:Ethnicity

####*********************************Significance Analysis of MocroArrays#################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
library(samr)  #
library(DT)
library(limma)
# source('http://bioconductor.org/biocLite.R') # Import biocLite() function into R environment biocLite('limma')
library(limma)
##1_control
##2_HAD
##3_HIV
##4_HIVE
SAMMultiClass_data_BG <- list(x =t(newpredictor_trainTestcombined_BG[,-c(711,712)]), y = t(targetAll), geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)

samr.objMultiClass_BG <- samr(SAMMultiClass_data_BG,resp.type="Multiclass",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.objMultiClass_BG)
delta.tableMultiClass_BG<-samr.compute.delta.table(samr.objMultiClass_BG)
delta.tableMultiClass_BG
del<-0.0594
siggenes.table_BG<- samr.compute.siggenes.table(samr.objMultiClass_BG,del,SAMMultiClass_data_BG,delta.tableMultiClass_BG,min.foldchange=2)
siggenes.table_BG
samr.plot(samr.objMultiClass_BG, del, min.foldchange=2)
pv=samr.pvalues.from.perms(samr.objMultiClass_BG$tt,samr.objMultiClass_BG$ttstar)
pv
write.table(siggenes.table_BG$genes.up,"DataFiles/BGMultiClassup.txt")
write.table(siggenes.table_BG$genes.lo,"DataFiles/BGMultiClasslow.txt")


#######Pairwise SAM Analysis######
##Control vs HAD
BG_y12 <- c(rep(1, 5), rep(2, 5),rep(1, 1), rep(2, 2))
##Control vs HIV
##Actual BG_y13 <- c(rep(1, 5), rep(3, 5),rep(1, 1), rep(3, 1))
BG_y13 <- c(rep(1, 5), rep(2, 5),rep(1, 1), rep(2, 1))
##Control vs HIVE
##Actual BG_y14 <- c(rep(1, 5), rep(4, 4),rep(1, 5), rep(4, 1))
BG_y14 <- c(rep(1, 5), rep(2, 4),rep(1, 1), rep(2, 1))


##HIV vs HAD
#Actual BG_y32 <- c(rep(3, 5), rep(2, 5),rep(3, 1), rep(2, 2))
BG_y32 <- c(rep(1, 5), rep(2, 5),rep(1, 1), rep(2, 2))
##HAD vs HIVE
##Actual BG_y24 <- c(rep(2, 5), rep(4, 4),rep(2, 2), rep(4, 1))
BG_y24 <- c(rep(1, 5), rep(2, 4),rep(1, 2), rep(2, 1))
##HIV vs HIVE
##Actual BG_y34 <- c(rep(3, 5), rep(4, 4),rep(3, 1), rep(4, 1))
BG_y34 <- c(rep(1, 5), rep(2, 4),rep(1, 1), rep(2, 1))

SAM_data_BG12 <- list(x =t(newpredictor_trainTestcombined_BG[c(1:5,11:15,20,22:23),-c(711,712)]), y = t(BG_y12), geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)
SAM_data_BG13 <- list(x =t(newpredictor_trainTestcombined_BG[c(1:5,6:10,20,21),-c(711,712)]), y = t(BG_y13),geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)
SAM_data_BG14 <- list(x =t(newpredictor_trainTestcombined_BG[c(1:5,16:19,20,24),-c(711,712)]), y = t(BG_y14),geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)

SAM_data_BG32 <- list(x =t(newpredictor_trainTestcombined_BG[c(6:10,11:15,21,22:23),-c(711,712)]), y = t(BG_y32),geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)
SAM_data_BG24 <- list(x =t(newpredictor_trainTestcombined_BG[c(11:15,16:19,22:23,24),-c(711,712)]), y = t(BG_y24),geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)
SAM_data_BG34 <- list(x =t(newpredictor_trainTestcombined_BG[c(6:10,16:19,21,24),-c(711,712)]), y = t(BG_y34),geneid = geneID_BG, genenames = genenames_BG, logged2 = TRUE)


##Checking for CIRBP geneID
which(genenames_BG=="BTN3A3") ##40
geneID[187]  ##1153

##1. Control vs HAD
samr.obj_BG12 <- samr(SAM_data_BG12, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_BG12)
delta.table_BG12 <- samr.compute.delta.table(samr.obj_BG12, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_BG12
del<-0.3295
siggenes.table_BG12<- samr.compute.siggenes.table(samr.obj_BG12,del,SAM_data_BG12,delta.table_BG12,min.foldchange=2)
siggenes.table_BG12
write.table(siggenes.table_BG12$genes.up,"DataFiles/BG12up.txt")
write.table(siggenes.table_BG12$genes.lo,"DataFiles/BG12low.txt")
samr.plot(samr.obj_BG12, del, min.foldchange=2)

##2. Control vs HIV
samr.obj_BG13 <- samr(SAM_data_BG13, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_BG13)
delta.table_BG13 <- samr.compute.delta.table(samr.obj_BG13, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_BG13
del<-0.4428
siggenes.table_BG13<- samr.compute.siggenes.table(samr.obj_BG13,del,SAM_data_BG13,delta.table_BG13,min.foldchange=2)
siggenes.table_BG13
write.table(siggenes.table_BG13$genes.up,"DataFiles/BG13up.txt")
write.table(siggenes.table_BG13$genes.lo,"DataFiles/BG13low.txt")
samr.plot(samr.obj_BG13, del, min.foldchange=2)

##3. Control vs HIVE
samr.obj_BG14 <- samr(SAM_data_BG14, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_BG14)
delta.table_BG14 <- samr.compute.delta.table(samr.obj_BG14, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_BG14
del<-0.7833
siggenes.table_BG14<- samr.compute.siggenes.table(samr.obj_BG14,del,SAM_data_BG14,delta.table_BG14,min.foldchange=2)
siggenes.table_BG14
write.table(siggenes.table_BG14$genes.up,"DataFiles/BG14up.txt")
write.table(siggenes.table_BG14$genes.lo,"DataFiles/BG14low.txt")
samr.plot(samr.obj_BG14, del, min.foldchange=2)

##4.HIV vs HAD
samr.obj_BG32 <- samr(SAM_data_BG32, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_BG32)
delta.table_BG32 <- samr.compute.delta.table(samr.obj_BG32, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_BG32
del<- 0.1000
siggenes.table_BG32<- samr.compute.siggenes.table(samr.obj_BG32,del,SAM_data_BG32,delta.table_BG32,min.foldchange=2)
siggenes.table_BG32
write.table(siggenes.table_BG32$genes.up,"DataFiles/BG32up.txt")
write.table(siggenes.table_BG32$genes.lo,"DataFiles/BG32low.txt")
samr.plot(samr.obj_BG32, del, min.foldchange=2)

##5.HAD vs HIVE
samr.obj_BG24 <- samr(SAM_data_BG24, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_BG24)
delta.table_BG24 <- samr.compute.delta.table(samr.obj_BG24, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_BG24
del<-0.2519
siggenes.table_BG24<- samr.compute.siggenes.table(samr.obj_BG24,del,SAM_data_BG24,delta.table_BG24,min.foldchange=2)
siggenes.table_BG24
write.table(siggenes.table_BG24$genes.up,"DataFiles/BG24up.txt")
write.table(siggenes.table_BG24$genes.lo,"DataFiles/BG24low.txt")
samr.plot(samr.obj_BG24, del, min.foldchange=2)

##6.HIV vs HIVE
samr.obj_BG34 <- samr(SAM_data_BG34, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_BG34)
delta.table_BG34 <- samr.compute.delta.table(samr.obj_BG34, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_BG34
del<-0.0979
siggenes.table_BG34<- samr.compute.siggenes.table(samr.obj_BG34,del,SAM_data_BG34,delta.table_BG34,min.foldchange=2)
siggenes.table_BG34
write.table(siggenes.table_BG34$genes.up,"DataFiles/BG34up.txt")
write.table(siggenes.table_BG34$genes.lo,"DataFiles/BG34low.txt")
samr.plot(samr.obj_BG34, del, min.foldchange=2)

save.image("HIVDataTrainTestDataCreation.RData")

#######################################################
load("HIVDataTrainTestDataCreation.RData")
###################
##
splom(t(newpredictor_trainTestcombinedAll), panel = panel.smoothScatter, raster = TRUE)

##Correlation Analysis for White MAtter
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)

hmcol <- colorRampPalette(brewer.pal(10,"PuOr"))(342)
newpredictor_trainTestcombined_BG[c(1:5,16:19,20,24),-c(711,712)]
corrs <- cor(t(newpredictor_trainTestcombined_BG[,-c(711,712)]), method = "spearman")
heatmap.2(corrs, margins = c(6, 11),col = hmcol, trace = "none", main="Correlation Analysis With Basal Ganglia",keysize = 1,Rowv=FALSE,Colv=FALSE,dendrogram='none', 
          srtCol = 56)


###Principal Component Analysis for Basal Ganglia
prcomps_BG<- prcomp(t(newpredictor_trainTestcombined_BG[,-c(711,712)]))
summary(prcomps_BG)
plot(prcomps_BG$rotation, type = "n",main="Principal Components for Basal Ganglia")
text(prcomps_BG$rotation, rownames(prcomps_BG$rotation))

install.packages("rgl")
library(rgl)
plot3d(prcomps_BG$rotation[,1:3],col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)))###loadings

install.packages("factoextra")
library(factoextra)
fviz_pca_ind(prcomps_BG,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
colors<-rainbow(4)
#colors <- c('#999999', '#E69F00', '#56B4E9','#1004F0')
colors <- colors[as.numeric(targetAll)]


###col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)) , Here the colors are assigned in the order the groups apppear in TargetAll data vector

plot3d(prcomps_BG$rotation[,1:3],col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)),type="s",image=TRUE,xlab="PC1",ylab="PC2",zlab="PC3")###loadings
legend3d("topright", legend = levels(targetAll), pch = 16, col = levels(targetAll), cex=1, inset=c(0.01))


snapshot3d(filename = "GraphsPAnalysis3DBG.png", fmt = 'png')


###Verimax Rotation on White Matter
rawLoadings_BG     <- prcomps_BG$rotation[,1:3] %*% diag(prcomps$sdev, 3, 3)
rotatedLoadings_BG<- varimax(rawLoadings_BG)$loadings
rotatedLoadings_BG

colors<-rainbow(4)
#colors <- c('#999999', '#E69F00', '#56B4E9','#1004F0')
colors <- colors[as.numeric(targetAll)]

plot3d(rotatedLoadings_BG[,1:3],col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)),type="s",image=TRUE,xlab="R1",ylab="R2",zlab="R3")###loadings
legend3d("topright", legend = levels(targetAll), pch = 20, col = levels(targetAll), cex=1, inset=c(0.01))

snapshot3d(filename ="RotatedPAnalysis3DBG.png", fmt = 'png')

########################################################################################
##Analysing Differentially Expressed Genes in Basal Ganglia with Moderated t-test
########################################################################################
design1 <- model.matrix(~0+targetAll)
colnames(design1)<-c("Control","HAD","HIV","HIVE")
array_BG<-arrayWeights(t(newpredictor_trainTestcombined_BG[,-c(711,712)]),design=design1)
fit1_BG <- lmFit(t(newpredictor_trainTestcombined_BG[,-c(711,712)]), design1,weights=array_BG)
dim(design1)
dim(fit1_BG)
contrast.matrix<-makeContrasts(HAD-Control,HIV-Control,HIVE-Control,HAD-HIV,HIVE-HAD,HIVE-HIV,levels=design1)

fit1_2_BG<-contrasts.fit(fit1_BG,contrast.matrix)

fit1_2_BG <- eBayes(fit1_2_BG)
BGT1<-topTable(fit1_2_BG,coef=1,adjust="BH",n="Inf") # A list of top genes differential expressed in HAD vs control
write.table(BGT1,"BGT1.txt")
BGT2<-topTable(fit1_2_BG,coef=2,adjust="BH",n="Inf") # A list of top genes differential expressed in HIV vs control
write.table(BGT2,"BGT2.txt")
BGT3<-topTable(fit1_2_BG,coef=3,adjust="BH",n="Inf") # A list of top genes differential expressed in HIVE vs control
write.table(BGT3,"BGT3.txt")
BGT4<-topTable(fit1_2_BG,coef=4,adjust="BH",n="Inf") # A list of top genes differential expressed in HAD vs HIV
write.table(BGT4,"BGT4.txt")
BGT5<-topTable(fit1_2_BG,coef=5,adjust="BH",n="Inf") # A list of top genes differential expressed in HIVE vs HAD
write.table(BGT5,"BGT5.txt")
BGT6<-topTable(fit1_2_BG,coef=6,adjust="BH",n="Inf") # A list of top genes differential expressed in HIVE vs HIV
write.table(BGT6,"BGT6.txt")


## we have used the cut-off in model after revision
results_BG = decideTests(fit1_2_BG, adj.P.Val=0.01)
write.table(results_BG,"BG_eBayes.txt")

##Venn Diagram Creation for Up and Down Genes
vennDiagram(results_BG[,c(1:3)],main="DE Genes in Basal Ganglia",include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))
#################################################################################

volcanoplot(fit1_2_BG)

#results = decideTests(fit2_2, adj.P.Val=0.05)
#vennDiagram(results)

plotMD(fit1_2_BG,coef=3,status =results1[,3],values=c(-1,1),hl.col=c("red","blue"))


#***********************************************************************************


##Checking for RBM3 geneID
which(genenames=="RBM3")# 1229
geneID[1229] ## 5935

##Checking for CIRBP geneID
which(genenames_BG=="CIRBP") ##187
geneID[187]  ##1153

##To identify the DE genes in all three sectors/all tw in BG
which(results3_1[,1]==1 & results3_1[,2]==1 & results3_1[,3]==1)
which(results3_1[,1]==1 & results3_1[,3]==1) ## HAD - Control and HIVE - Control


save.image("BasalGangliaAnalysis.RData")
