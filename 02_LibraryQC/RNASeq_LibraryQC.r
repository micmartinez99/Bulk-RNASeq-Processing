
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


##EDIT LINE BELOW:  SET this to where the output is desired
files <- list.files("/Users/mikemartinez/Desktop/WT_vs_Control/Counts/", pattern = "*.tsv$")

# Check that all files are being listed as expected
files

# Initialize a counter
i<-0

# Iteratively read in the files
while (i<length(files)){
  if (i == 0){
    print(paste0("Read in file : ",files[i+1]))
    fileName=paste0(countsDir,"/",files[i+1])
    
    # Read in the first file
    df1<-read.table(fileName,sep="\t",stringsAsFactors = FALSE,header=F)
    colnames(df1)<-c("geneID",strsplit(files[i+1],".tsv"))
    
    # Increment the counter
    i<-i+1
  }else{
  
  	# Read in the rest of the files
    print(paste0("Read in file : ",files[i+1]))
    fileName=paste0(countsDir,"/",files[i+1])
    df2<-read.table(fileName,sep="\t",stringsAsFactors = FALSE,header=F)
    colnames(df2)<-c("geneID",strsplit(files[i+1],".tsv"))
    
    # Concatenate the counts file as it is read in
    df1<-merge(df1,df2,by.x="geneID",by.y="geneID",sort=FALSE)
    i<-i+1
  }  
}

# Check that the file looks okay
head(df1)

# Check that ncol == number of samples you have
dim(df1)
df1<-df1[(!stringr::str_starts(df1[["geneID"]],"__")),]
dim(df1)

df1<-column_to_rownames(df1,var = "geneID")
head(df1)

# Pre-filter low counts
df1<-df1[(rowSums(df1)>10),]

# Save file and adjust as needed
write.csv(df1, file = "WT_vs_ControlNormal_Counts.csv")


# Library QC
df.m <- reshape2::melt(df1, id.vars =NULL)

# Plot to asses slibrary sizes
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

