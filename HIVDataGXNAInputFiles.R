##This is the File 4 for All three Analysis WM, BG and FC

##Creating GXNA Input Files
setwd("MD/HIVData/ModelAfterRevision")

load("ModeAfterRevisionGXNAInputCreation.RData")
load("HIVDataTrainTestDataCreation.RData")

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GEOquery")
BiocManager::install("affy", version = "3.8")
BiocManager::install(c("gcrma","org.Hs.eg.db"))
BiocManager::install("hgu133plus2.db")

##To run the following package, First Download the tar file from 
###https://bioconductor.org/packages/release/bioc/html/BiocGenerics.html
##And then install the package locally.
BiocManager::install("BiocGenerics")
###First Install the following packages using Install packages from local files: from "G:\MD"" 
##original source files are available at http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp

BiocManager::install(c("hgu133plus2hsentrezgcdf","hgu133plus2hsentrezgprobe","hgu133plus2hsentrezg.db"))
#BiocManager::install(c("hgu133ahsentrezgcdf","hgu133ahsentrezgprobe","hgu133ahsentrezg.db"))

library(GEOquery)
library(affy)
library(gcrma)
library(AnnotationDbi)
library(hgu133plus2hsentrezgcdf) #cdfname="HGU133Plus2_Hs_ENTREZG"
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezg.db)

####White Matter GXNA File Creations

trimNumber = function(x) {
  return(sprintf("%.3f", x)) # keep only three decimals
}

writeArray = function(arrayName, fileName) {
  write.table(arrayName, file = fileName,
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}
##Making a prbname_WM as a charactor variable 
prbname_WM<-as.character(prbname_WM)
writeExpression = function(MAListName, fileName) {
  writeArray(cbind(prbname_WM, 
                   apply(filt_TrainData_WM, 2, trimNumber)), fileName)
}

prepareGXNA = function(MAListName, phenotype, projectName) {
  expFile = paste(projectName, "exp", sep = ".")
  phenoFile = paste(projectName, "phe", sep = ".")
  writeExpression(MAListName, expFile)
  writeArray(phenotype, phenoFile) # must transpose phenotype array
}

prepareGXNA(filt_TrainData_WM, trainY, "HIV_WM_AfterRevision")

##In case, this also works
writeArray(cbind(prbname_WM, 
                 apply(filt_TrainData_WM, 2, trimNumber)), "HIV_WM_AfterRevision.exp")

###This Worked to create the annotion file#####
##This works for all three BG, WM and FC
##Hence,  need to be created only once.

library("hgu133plus2hsentrezgcdf", character.only=TRUE) # load package
map = as.list(get(paste("hgu133plus2hsentrezg","ENTREZID",sep = "")))
  mapFile = paste("hgu133plus2hsentrezg", "ann", sep = ".")
  write.table(t(rbind(names(map), map)), file = mapFile,
              row.names=FALSE, col.names=FALSE)

#########################################  
####Basal Ganglia GXNA File Creations
########################################  
  trimNumber = function(x) {
    return(sprintf("%.3f", x)) # keep only three decimals
  }
  
  writeArray = function(arrayName, fileName) {
    write.table(arrayName, file = fileName,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  ##Making a prbname_WM as a charactor variable 
  prbname_BG<-as.character(prbname_BG)
  writeExpression = function(MAListName, fileName) {
    writeArray(cbind(prbname_BG, 
                     apply(filt_TrainData_BG, 2, trimNumber)), fileName)
  }
  
  prepareGXNA = function(MAListName, phenotype, projectName) {
    expFile = paste(projectName, "exp", sep = ".")
    phenoFile = paste(projectName, "phe", sep = ".")
    writeExpression(MAListName, expFile)
    writeArray(phenotype, phenoFile) # must transpose phenotype array
  }
  
  prepareGXNA(filt_TrainData_BG, trainY, "HIV_BG_AfterRevision")
  
  ##In case, this also works
  writeArray(cbind(prbname_BG, 
                   apply(filt_TrainData_BG, 2, trimNumber)), "HIV_BG_AfterRevision.exp")
  


  #########################################  
  ####Frontal Cortex GXNA File Creations
  ########################################  
  trimNumber = function(x) {
    return(sprintf("%.3f", x)) # keep only three decimals
  }
  
  writeArray = function(arrayName, fileName) {
    write.table(arrayName, file = fileName,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  ##Making a prbname_WM as a charactor variable 
  prbname_FC<-as.character(prbname_FC)
  writeExpression = function(MAListName, fileName) {
    writeArray(cbind(prbname_FC, 
                     apply(filt_TrainData_FC, 2, trimNumber)), fileName)
  }
  
  prepareGXNA = function(MAListName, phenotype, projectName) {
    expFile = paste(projectName, "exp", sep = ".")
    phenoFile = paste(projectName, "phe", sep = ".")
    writeExpression(MAListName, expFile)
    writeArray(phenotype, phenoFile) # must transpose phenotype array
  }
  
  prepareGXNA(filt_TrainData_FC, trainY, "HIV_FC_AfterRevision")
  
  ##In case, this also works
  writeArray(cbind(prbname_FC, 
                   apply(filt_TrainData_FC, 2, trimNumber)), "HIV_BG_AfterRevision.exp") 
 
save.image("ModeAfterRevisionGXNAInputCreation.RData")

###############################################

###Part II GXNA Analysis#############
#########################################
#mapNewData<-read.table("prbEntrezId.txt",header=TRUE,sep="\t")
#mapNewData <- as.data.frame(sapply(mapNewData, function(x) gsub("\"", "", x)))
#   mapFile = paste("hgu133plus2", "ann", sep = ".")
#  write.table(mapNewData,file = mapFile,row.names=FALSE, col.names=FALSE,quote=FALSE,eol = "")


#To Run GXNA software
#Open command prompt
#Type "cd Desktop"
#cd gxna
#C:\Users\18134\Desktop>cd gxna


###Do the file creation .phe,.ann. and .exp file creations based on the GXNA tool 
###Save the files at C:\Users\18134\Desktop\gxna
##load the mm vector from the R file called "HIVDataGXNAInputFilesR.R" saved at C:\Users\18134\Desktop\gxna

### Load the initial clusters from GXNA into an array called "mm_WM[[]]"
load("ModeAfterRevisionGXNAInputCreation.RData")

##Finally we run the RF algorithm. NOTE: we use an ODD number for ntree. This is because when the forest/ensembl is used on test data, ties are broken randomly. 
##Having an odd number of trees avoids this issue and makes the model fully deterministic. 
##Also note, we will use down-sampling to attempt to compensate for unequal class-sizes 
#(less ##relapses than non-relapses).

###########################
##White Matter
##########################

##Path Before Model After Revision
path = "C:/Users/18134/Desktop/gxna/Files/results"

##Path After Model Revision
path = "C:/Users/18134/Desktop/gxna"
out.file<-""
file.names <- dir(path, pattern =".txt")
for(i in 1:3){
  file <- read.table(file.names[i],header=TRUE, sep=";", stringsAsFactors=FALSE)
  out.file <- rbind(out.file, file)
  }


###Run the Following Code to import the data into mm matrix. This is for Model After Revision
##Data from HIV_WM_AfterRevision
setwd("C:/Users/18134/Desktop/gxna")
ninitClusters=50
ind=ninitClusters-1
mm_WM<-matrix(list(), 1,ninitClusters)
for(i in 1:ninitClusters){
mapfile=paste("HIV_WM_AfterRevision_1_",i-1,".txt",sep="")
mm_WM[[i]]=read.table(mapfile)}

save.image("WM_gxnaTrainClusters.RData")

###For HIVFull Data 
setwd("C:/Users/18134/Desktop/gxna")
ninitClusters=50
ind=ninitClusters-1
mm<-matrix(list(), 1,ninitClusters)
for(i in 1:ninitClusters){
  mapfile=paste("HIVFull_001_",i-1,".txt",sep="")
  mm[[i]]=read.table(mapfile)}

###For HIVFull Data 
setwd("C:/Users/18134/Desktop/gxna")
ninitClusters=50
ind=ninitClusters-1
mm2<-matrix(list(), 1,ninitClusters)
for(i in 1:ninitClusters){
  mapfile=paste("HIVFull_002_",i-1,".txt",sep="")
  mm2[[i]]=read.table(mapfile)}
save.image("gxnaTrainClusters.RData")


###########################
##Basal Ganglia
##########################

##Path Before Model After Revision
path = "C:/Users/18134/Desktop/gxna/Files/results"


###Run the Following Code to import the data into mm_BG matrix. This is for Model After Revision
##Data from HIV_BG_AfterRevision
setwd("C:/Users/18134/Desktop/gxna")
ninitClusters=50
ind=ninitClusters-1
mm_BG<-matrix(list(), 1,ninitClusters)
for(i in 1:ninitClusters){
  mapfile=paste("HIV_BG_AfterRevision_1_",i-1,".txt",sep="")
  mm_BG[[i]]=read.table(mapfile)}

save.image("BG_gxnaTrainClusters.RData")


###########################
##Frontal Cortex
##########################

##Path Before Model After Revision
path = "C:/Users/18134/Desktop/gxna/Files/results"


###Run the Following Code to import the data into mm_FC matrix. This is for Model After Revision
##Data from HIV_FC_AfterRevision
setwd("C:/Users/18134/Desktop/gxna")

ninitClusters=50
ind=ninitClusters-1
mm_FC<-matrix(list(), 1,ninitClusters)
for(i in 1:ninitClusters){
  mapfile=paste("HIV_FC_AfterRevision_1_",i-1,".txt",sep="")
  mm_FC[[i]]=read.table(mapfile)}

save.image("FC_gxnaTrainClusters.RData")

