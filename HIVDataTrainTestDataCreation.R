##DataFile: file:///C:/Users/18134/OneDrive%20-%20The%20University%20of%20Texas-Rio%20Grande%20Valley/data2/GEO%20DataSet%20Browser.html
###To get custom  CDf file: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
##ftp://ftp.ncbi.nlm.nih.gov/geo/
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7146 ##Diabetes Data
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35864 ##HIV Data
##https://nntc.org/content/gene_array
###SAM Analysis: https://mdozmorov.github.io/BIOS567/assets/presentation_diffexpression/DiffExpr_SAM.html
##Gene Filtering http://math.usu.edu/jrstevens/stat5570/3.2.Filtering.pdf

##This is the HIV File 1

#Part I
##HIV File 1 - Part II is down in this same file

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())


chooseCRANmirror()

tempdir="G:/testData"
outdir="G:/MD/HIVData"

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GEOquery")
BiocManager::install("beadarray")
BiocManager::install("affy", version = "3.8")
BiocManager::install(c("gcrma","org.Hs.eg.db"))
BiocManager::install("hgu133plus2.db")

###First Install the following packages using Install packages from local files: from "G:\MD"" 
##original source files are available at http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp

BiocManager::install(c("hgu133plus2hsentrezgcdf","hgu133plus2hsentrezgprobe","hgu133plus2hsentrezg.db"))
#BiocManager::install(c("hgu133ahsentrezgcdf","hgu133ahsentrezgprobe","hgu133ahsentrezg.db"))

library(GEOquery)
library(affy)
library(gcrma)
library(hgu133plus2hsentrezgcdf) #cdfname="HGU133Plus2_Hs_ENTREZG"
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezg.db)

#Set a data directory, download a GEO dataset, unpack and gunzip, 
#and create a list of files for processing.


####Creating the train data set
setwd(tempdir)
getGEOSuppFiles("GSE35864")

setwd(paste(tempdir,"GSE35864", sep="/"))
untar("GSE35864_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

#Create an AffyBatch object and perform GCRMA normalization.
setwd(paste(tempdir,"GSE35864/data", sep="/"))
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="HGU133Plus2_Hs_ENTREZG") 



data.gcrma.norm=gcrma(raw.data)

##look for probe names starting with "AFFX". These are the control genes.
# dim(gcrmaOLD)=20481 by 72
gcrmaOLD=exprs(data.gcrma.norm)


# dim(gcrma)=20419 by 72
gcrma=gcrmaOLD[which(!grepl("AFFX", rownames(gcrmaOLD))),]
#gcrmacontrol=gcrmaOLD[20420:20481,]

#Map probe set identifiers to Entrez gene symbols and IDs and then combine with raw data.
#ls("package:hgu133plus2hsentrezg.db") 
probes=row.names(gcrma)
symbol = unlist(mget(probes, hgu133plus2hsentrezgSYMBOL))
ID = unlist(mget(probes, hgu133plus2hsentrezgENTREZID))
gcrma=cbind(probes,ID,symbol,gcrma)ls("package:hugene10stprobeset.db")
select(hgu133plus2hsentrezg.db,c("10006_at","100128079_at","100128537_at"),c("SYMBOL","ENTREZID","GENENAME"))
#Write GCRMA-normalized and mapped data to file.
setwd(outdir)
write.table(gcrma, file = "HIVFullData_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Get the clinical details for this dataset.

gse <- getGEO("GSE35864", GSEMatrix = TRUE)
show(gse)
dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])

GSE35864_clindata=pData(gse[[1]])
write.table(GSE35864_clindata, "HIVFullDataset_clindetails.txt", quote = FALSE, sep = "\t", row.names = FALSE)
##Just uploading the only seleceted variable information from the clinical datafile
GSE35864_clindata2=pData(gse[[1]])[,c(2,6,37:41)]
write.table(GSE35864_clindata2, "HIVFullDataset_clindetails2.txt", quote = FALSE, sep = "\t", row.names = FALSE)

########################################################################

## HIV File Analysis Number 1 - Part II


# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

#### Creating the Random Forest Model###

install.packages("ROCR")
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("checkmate")
BiocManager::install("genefilter")


library(ROCR)
library(genefilter)
library(ggplot2)
library(checkmate)
library(Hmisc)


datafile="HIVFullData_gcrma.txt" 
clindatafile="HIVFullDataset_clindetails2.txt"

outfile="trainset_RFoutput.txt"
varimp_pdffile="trainset_varImps.pdf"
MDS_pdffile="trainset_MDS.pdf"
ROC_pdffile="trainset_ROC.pdf"
case_pred_outfile="trainset_CasePredictions.txt"
vote_dist_pdffile="trainset_vote_dist.pdf"

data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t") # dim(data_import)=20419 by 75
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

clin_data_order=order(clin_data_import[,"geo_accession"])
clindata=clin_data_import[clin_data_order,]

#Order data without first three columns,then add 3 to get correct index in original file
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 

#grab first three columns, and then remaining columns in order determined above
rawdata=data_import[,c(1:3,data_order)] #dim(rawdata)= 20419  by 75
header=colnames(rawdata)
rawdata=rawdata[which(!is.na(rawdata[,3])),] #Remove rows with missing gene symbol; dim(rawdata)= 20408  by  75

##We will assign the variables to a new data structure. 
##Extract just the expression values from the filtered data and transpose the matrix. 
##The latter is necessary because RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. 
##Finally, assign gene symbol as the predictor name.

#predictor_data=t(rawdata[,4:length(header)])
predictor_data=t(rawdata[,4:length(header)])


predictor_names=c(as.vector(rawdata[,3])) #gene symbol
predictor_probenames=c(as.vector(rawdata[,1])) #Probe ID
#colnames(predictor_data)=predictor_names
colnames(predictor_data)=predictor_probenames

#Get potential predictor variables for brain tissue types "White Matter","Frontal Cortex" and "Basal Ganglia".
predictor_dataWM=t(rawdata[,4:27])
predictor_dataFC=t(rawdata[,28:51])
predictor_dataBG=t(rawdata[,52:75])

GXNAFile1=rawdata[,1:3]
##As a final step before the Random Forest classification, we have to set the variable we are trying to predict as our ##target variable. In this case, it is relapse status.

target=c(rep("control",6),rep("HIV",6),rep("HAD",7),rep("HIVE",5))
target=as.factor(target)

#target=c(rep("control",6),rep("HIV",6),rep("HAD",7),rep("HIVE",5),rep("control",6),rep("HIV",6),rep("HAD",7),rep("HIVE",5),rep("control",6),rep("HIV",6),rep("HAD",7),rep("HIVE",5))
#target=as.factor(target)


####Creating the training and Testing datasets
install.packages("caret")
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(target, p = .70, list = FALSE, times = 1)
head(trainIndex)

Train_WM <- predictor_dataWM[ trainIndex,]
Test_WM <- predictor_dataWM[-trainIndex,]

Train_FC <- predictor_dataFC[trainIndex,]
Test_FC <- predictor_dataFC[-trainIndex,]

Train_BG <- predictor_dataBG[ trainIndex,]
Test_BG <- predictor_dataBG[-trainIndex,]

dim(Train_WM)
dim(Test_FC)

library(randomForest) 
set.seed(1234) 

####No need the following####
trainX_WM=t(Train_WM)
testX_WM=t(Test_WM)

trainX_FC=t(Train_FC)
testX_FC=t(Test_FC)

trainX_BG=t(Train_BG)
testX_BG=t(Test_BG)

##Only one trainY is needed as they are same for all three groups
trainY=target[trainIndex]  #dim(trainX)=20408 by 19
testY=target[-trainIndex]

#####Taking the Demographic variables Age, Race and Tissue type

clinicalData<-read.table("HIVFullDataset_clindetails2.txt",header = TRUE,sep="\t")
clincData24Individual<-clinicalData[1:24,]
TrainXDem=clincData24Individual[trainIndex,c(3,6,7)] ##dim(TrainXDem)= 19 by  3
TestXDem=clincData24Individual[-trainIndex,c(3,6,7)] ##dim(TestXDem)= 5 by  3

##*****************
##GCRMA normlisatiobn return expression value in log2 scale and by default will perfromed the quantile normalization

source("https://bioconductor.org/biocLite.R")
biocLite("geneplotter")
biocLite("RColorBrewer")
# make plots
library(geneplotter)
library(RColorBrewer)
blues.ramp<-colorRampPalette(brewer.pal(9,"Blues")[-1])
dCol<-densCols(log(gene.mean),log(gene.sd),colramp=blues.ramp)

##Take the un-logged data to see the effect of the filter on White Matter
unloggedTrainX_WM=2^trainX_WM
gene.mean_WM<-apply(unloggedTrainX_WM,1,mean)
gene.sd_WM <-apply(unloggedTrainX_WM,1,sd)
gene.cv_WM <-gene.sd_WM/gene.mean_WM

dev.off()
par(mfrow=c(1,2))
plot(gene.mean_WM,gene.sd_WM,log='xy',col=dCol,pch=16,cex=0.1,main="Look at Mean White Matter")
abline(v=100,lwd=3,col='red')

hist(log2(gene.cv_WM),main="Look at CV White Matter")
abline(v=log(.7),lwd=3,col='red')

##Take the un-logged data to see the effect of the filter for Frontal Cortex
unloggedTrainX_FC=2^trainX_FC
gene.mean_FC<-apply(unloggedTrainX_FC,1,mean)
gene.sd_FC <-apply(unloggedTrainX_FC,1,sd)
gene.cv_FC <-gene.sd_FC/gene.mean_FC

dev.off()
par(mfrow=c(1,2))
plot(gene.mean_FC,gene.sd_FC,log='xy',col=dCol,pch=16,cex=0.1,main="Look at Mean Frontal Cortex")
abline(v=100,lwd=3,col='red')

hist(log2(gene.cv_FC),main="Look at CV Frontal Cortex")
abline(v=log(.7),lwd=3,col='red')

dev.off()
hist(unloggedTrainX_FC,xlim=c(0,10000))
dev.off()
hist(trainX_FC)
#> max(trainX_FC)
#15.46848
# min(trainX_FC)
#2.170877
##Take the un-logged data to see the effect of the filter for Basal Ganglia
unloggedTrainX_BG=2^trainX_BG
hist(unloggedTrainX_BG, col = "gray", main="Basal Ganglia")
gene.mean_BG<-apply(unloggedTrainX_BG,1,mean)
gene.sd_BG <-apply(unloggedTrainX_BG,1,sd)
gene.cv_BG <-gene.sd_BG/gene.mean_BG

dev.off()
par(mfrow=c(1,2))
plot(gene.mean_BG,gene.sd_BG,log='xy',col=dCol,pch=16,cex=0.1,main="Look at Mean of Basal Ganglia")
abline(v=100,lwd=3,col='red')

hist(log2(gene.cv_BG),main="Look at CV of Basal Ganglia")
abline(v=log(.7),lwd=3,col='red')


###### All Plots in one Graph

dev.off()
par(mfrow=c(3,2))
plot(gene.mean_WM,gene.sd_WM,log='xy',col=dCol,pch=16,cex=0.1,main="Look at Mean White Matter")
abline(v=100,lwd=3,col='red')

hist(log2(gene.cv_WM),main="Look at CV White Matter")
abline(v=log2(.7),lwd=3,col='red')

plot(gene.mean_FC,gene.sd_FC,log='xy',col=dCol,pch=16,cex=0.1,main="Look at Mean Frontal Cortex")
abline(v=100,lwd=3,col='red')

hist(log2(gene.cv_FC),main="Look at CV Frontal Cortex")
abline(v=log2(.7),lwd=3,col='red')
plot(gene.mean_BG,gene.sd_BG,log='xy',col=dCol,pch=16,cex=0.1,main="Look at Mean of Basal Ganglia")
abline(v=100,lwd=3,col='red')

hist(log2(gene.cv_BG),main="Look at CV of Basal Ganglia")
abline(v=log2(.7),lwd=3,col='red')

###############################################################################
###Next we filter out any variables (genes) that are not expressed or 
##do not have enough variance to be informative in classification. 
##We will first take the values and un-log2 them, then filter out any genes according to following criteria 
#(recommended in multtest/MTP documentation): 
##(1) At least 20% of samples should have raw intensity greater than 100; 
##(2) The coefficient of variation (sd/mean) is between 0.7 and 10.


#### Creating the Random Forest Model###

install.packages("ROCR")
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("checkmate")
BiocManager::install("genefilter")


library(ROCR)
library(genefilter)
library(ggplot2)
library(checkmate)
library(Hmisc)

ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))


##Filtering genes for White Matter
filt_WM=genefilter(2^trainX_WM,ffun)
filt_TrainData_WM=trainX_WM[filt_WM,]  # dim(filt_TrainData_WM)= 1481 by  19
filt_TestData_WM=testX_WM[filt_WM,]  # dim(filt_TestData_WM)= 1481 by  5
filt_genes_WM=GXNAFile1[filt_WM,]
filteredDataAll_WM=predictor_dataWM[,filt_WM]#dim(filteredDataAll_WM)=24 by 1481
prbname_WM=rawdata[filt_WM,1]

##Filtering genes for Frontal Cortex
filt_FC=genefilter(2^trainX_FC,ffun)
filt_TrainData_FC=trainX_FC[filt_FC,]  # dim(filt_TrainData_FC)= 794 by  19
filt_TestData_FC=testX_FC[filt_FC,]  # dim(filt_TestData_FC)= 794 by  5
filt_genes_FC=GXNAFile1[filt_FC,]
filteredDataAll_FC=predictor_dataFC[,filt_FC]#dim(filteredDataAll_FC)=24 by 794
prbname_FC=rawdata[filt_FC,1]


##Filtering genes for Basal Ganglia
filt_BG=genefilter(2^trainX_BG,ffun)
filt_TrainData_BG=trainX_BG[filt_BG,]  # dim(filt_TrainData_BG)= 710 by  19
filt_TestData_BG=testX_BG[filt_BG,]  # dim(filt_TestData_BG)= 710 by  5
filt_genes_BG=GXNAFile1[filt_BG,]
filteredDataAll_BG=predictor_dataBG[,filt_BG]#dim(filteredDataAll_BG)=24 by 710
prbname_BG=rawdata[filt_BG,1]

save.image("HIVDataTrainTestDataCreation.RData")
###########################################################
load("HIVDataTrainTestDataCreation.RData")
###########################################################


###Just to see how the unloged filtered data looks for White Matter
filt_TrainDataunloged_WM=2^filt_TrainData_WM
gene.mean1_WM<-apply(filt_TrainDataunloged_WM,1,mean)
gene.sd1_WM <-apply(filt_TrainDataunloged_WM,1,sd)
gene.cv1_WM <-gene.sd1_WM/gene.mean1_WM
dev.off()
par(mfrow=c(1,2))
plot(gene.mean1_WM,gene.sd1_WM,log='xy',col=dCol,pch=16,cex=0.1,main="After Filtering expression for WM>100")
hist(log2(gene.cv1_WM),main="After Filtering CV 0.7-10 for WM")


###Just to see how the unloged filtered data looks for Frontal Cortex
filt_TrainDataunloged_FC=2^filt_TrainData_FC
gene.mean1_FC<-apply(filt_TrainDataunloged_FC,1,mean)
gene.sd1_FC <-apply(filt_TrainDataunloged_FC,1,sd)
gene.cv1_FC <-gene.sd1_FC/gene.mean1_FC

dev.off()
par(mfrow=c(1,2))
plot(gene.mean1_FC,gene.sd1_FC,log='xy',col=dCol,pch=16,cex=0.1,main="After Filtering expression for FC>100")
hist(log2(gene.cv1_FC),main="After Filtering CV 0.7-10 for FC")


###Just to see how the unloged filtered data looks for Basal Ganglia
filt_TrainDataunloged_BG=2^filt_TrainData_BG
gene.mean1_BG<-apply(filt_TrainDataunloged_BG,1,mean)
gene.sd1_BG <-apply(filt_TrainDataunloged_BG,1,sd)
gene.cv1_BG <-gene.sd1_BG/gene.mean1_BG

dev.off()
par(mfrow=c(1,2))
plot(gene.mean1_BG,gene.sd1_BG,log='xy',col=dCol,pch=16,cex=0.1,main="After Filtering expression for BG>100")
hist(log2(gene.cv1_BG),main="After Filtering CV 0.7-10 for BG")


##Three Unlogged plots for Data After Filtering 

dev.off()
par(mfrow=c(3,2))

plot(gene.mean1_WM,gene.sd1_WM,log='xy',col=dCol,pch=16,cex=0.1,main="After Filtering expression for WM>100")
hist(log2(gene.cv1_WM),main="After Filtering CV 0.7-10 for WM")

plot(gene.mean1_FC,gene.sd1_FC,log='xy',col=dCol,pch=16,cex=0.1,main="After Filtering expression for FC>100")
hist(log2(gene.cv1_FC),main="After Filtering CV 0.7-10 for FC")

plot(gene.mean1_BG,gene.sd1_BG,log='xy',col=dCol,pch=16,cex=0.1,main="After Filtering expression for BG>100")
hist(log2(gene.cv1_BG),main="After Filtering CV 0.7-10 for BG")


#############################################################################
save.image("HIVDataTrainTestDataCreation.RData")
################################################################################


