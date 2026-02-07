## This scripts is to perform topic analysis using pycisTopic for robust meta-regulatory programs
## The aim is to extract selected cells per class, conver to cistopic object,
## run topic analysis for a number of topics, then identify top region per topic
## splititng the data into two to obtain robust class specific topics
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/conda/envs/scenicplus/bin/python3.8")
use_python("/home/loc/to/conda/envs/scenicplus/bin/python3.8")
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(scater)
  
})

addArchRGenome("hg38")
addArchRThreads(threads = 7)
cls = commandArgs(trailingOnly=TRUE)
Split <- "B" # A or B
#directories
DIR_ATAC <- "~/ATACana/Outs/"
setwd(DIR_ATAC)
print(paste0("Performing cisTopic analysis for class: ",cls))
####### Processing ArchR object to obtain the peak matrix for selected cells ####### 
## 1. load metadata and select cells---------------
load("~/extrafiles/plotdataHB_ATACv6.RData")
pd <- plot.data.ATAC[plot.data.ATAC$Class==cls & plot.data.ATAC$SplitTopic==Split,]
rm(plot.data.ATAC)
#subset further, looks like 3k cells should be good
cells <- c()
cls <- unique(pd$Cluster_step2)
for(x in cls){
  cells_sel <- row.names(pd)[pd$Cluster_step2==x]
  if(length(cells_sel)>3000){
    cells_sel <- sample(cells_sel,3000)
    cells <- c(cells,cells_sel)
  } else {
    cells <- c(cells,cells_sel)
  }
  rm(cells_sel)
}
rm(x)

pd <- pd[cells,]
## further remove cells that are less than 10 cells
tb <- data.frame(table(pd$Cluster_step2))
tb <- as.character(tb$Var1)[tb$Freq>9]
pd <- pd[pd$Annotationlv2_step2 %in% tb,]
cells <- row.names(pd)

## 2. load ArchR project to obtain peak matrix-------------
proj <- loadArchRProject(path = "comATACx/")
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)
#robust peaks to subset to
peaks <- read.delim("PeakCom/AllPeaks_robust.fil.sorted.bed", header = FALSE)
peaks <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)

#process peak Matrix
PM <- assay(PeakMatrix)
region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)
row.names(PM) <- row.names(region_names)
PM[1:4,1:4]
PM <- PM[peaks,cells] ##important to check cell order
dim(PM)
## fix cell ids
row.names(pd) <- colnames(PM) <- pd$CellID

rm(PeakMatrix,peaks,region_names)

## 3. Processing inpyut for cisTopic conversion-------------
##removing peaks not present
rs <- rowSums(PM)
keep <- rs>10 ##filter peaks not present enough time
table(keep)
PM <- PM[keep,]
dim(PM)
peaks <- row.names(PM)
# adjust rownames of metadata and colnames of count matrix
a <- rownames(pd)
b <- paste(a,"___cisTopic",sep = "")
rownames(pd) <- b
head(pd)
out_loc <- paste0(args,"_",Split)
## 4. use python to generate cisTopic object ----------------
repl_python()

import os
import pandas as pd
import numpy
import warnings
import pycisTopic
from scipy import io
import pickle
from pycisTopic.cistopic_class import *

Clss = r.out_loc
projDir = os.path.join("TopicMP/",Clss)
if not os.path.exists(projDir):
  os.makedirs(projDir)

#
count_mat = sparse.csr_matrix(r.PM, dtype="int")
cn = r.a
peaks=r.peaks
#Create cisTopic object
path_to_blacklist='/home/ann/hg38-blacklist.v2_noncanchr.bed'
print("Creating pycistopic object from peak count matrix  ")

##path to fragment file not used
cistopic_obj = create_cistopic_object(fragment_matrix=count_mat,
                                      cell_names=cn, region_names= peaks,
                                      path_to_blacklist=path_to_blacklist)


print(" Adding metadata to cistopic object  ")

cistopic_obj.add_cell_data(r.pd)
print(cistopic_obj)
#Save
print(" Saving ...  ")
pickle.dump(cistopic_obj, open(os.path.join(projDir,Clss+'_cistopic_obj.pkl'), 'wb'))


exit
print(paste0("Finished creating cisTopic object for: ",out_loc))
q()
