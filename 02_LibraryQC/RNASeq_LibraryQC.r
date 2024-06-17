
rm(list=ls())
gc()

library(dplyr)
library(tidyverse)
library(DESeq2)
library("genefilter")
library("ggplot2")
library("grDevices")
library("ggrepel")
library("pheatmap")
library("RColorBrewer")
library("gplots")

#Set working directory
setwd <- "/Users/mikemartinez/Desktop/WT_vs_Control/Counts/"
countsDir <- "/Users/mikemartinez/Desktop/WT_vs_Control/Counts/"

# ref<-read.csv("/Users/astral/Documents/DataAnalysis/geneAnnotnMaster/MusMusculus_gtf_ensembl_infoMerged_reference.csv",
#               sep=",", header = T)
# head(ref)
# ref<-ref %>% column_to_rownames(var="geneID")
# head(ref)

##EDIT LINE BELOW:  SET this to where the output is desired
files <- list.files("/Users/mikemartinez/Desktop/WT_vs_Control/Counts/", pattern = "*.tsv$")
files
i<-0
while (i<length(files)){
  if (i == 0){
    print(paste0("Read in file : ",files[i+1]))
    fileName=paste0(countsDir,"/",files[i+1])
    df1<-read.table(fileName,sep="\t",stringsAsFactors = FALSE,header=F)
    colnames(df1)<-c("geneID",strsplit(files[i+1],".tsv"))
    i<-i+1
  }else{
    print(paste0("Read in file : ",files[i+1]))
    fileName=paste0(countsDir,"/",files[i+1])
    df2<-read.table(fileName,sep="\t",stringsAsFactors = FALSE,header=F)
    colnames(df2)<-c("geneID",strsplit(files[i+1],".tsv"))
    df1<-merge(df1,df2,by.x="geneID",by.y="geneID",sort=FALSE)
    i<-i+1
  }  
}
head(df1)
dim(df1)
df1<-df1[(!stringr::str_starts(df1[["geneID"]],"__")),]
dim(df1)

df1<-column_to_rownames(df1,var = "geneID")
head(df1)


df1<-df1[(rowSums(df1)>10),]
write.csv(df1, file = "WT_vs_ControlNormal_Counts.csv")


######COunts QC
df.m <- reshape2::melt(df1, id.vars =NULL)
QC <- ggplot(df.m,aes(factor(variable),log10(value),fill=variable)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  labs(title = "Library QC")
QC



samples<-colnames(df1)
sampleTableC<-read.csv("sampleTable.csv")
head(sampleTableC)

rownames(sampleTableC)<-sampleTableC$samples
sampleTableC

sum(samples %in% sampleTableC$samples)
if (sum(samples %in% sampleTableC$samples) == length(sampleTableC$samples)){
  message("Good News !!! Samples in count matrix matches with that of in sampleTable")
}

save(df1,sampleTableC,ref,file="input.rdata")

