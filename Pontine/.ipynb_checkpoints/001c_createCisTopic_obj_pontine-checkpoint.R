## This script converts ATAC data from ArchR into cisTopic format
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/path/to/conda/envs/scenicplus/bin/python3.8")
use_python("/path/to/conda/envs/scenicplus/bin/python3.8")
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(scater)
})
addArchRGenome("hg38")

DIR <- "~/ATACana/Outs/"
DIR_SCENIC <- "~/SCENICOut/PN/"

setwd(DIR)
##-----------
#Load the ATAC project
load("extrafiles/plotdataHB_ATACv6_0404255step2.RData")
pd <- data.frame(plot.data.ATAC)

cl <- c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES",
        "Neurons_I:PonNuc_def","Neurons_I:PonNuc_diff_early","Neurons_I:PonNuc_diff_late",
        "Neurons_II:PonNuc_mat")
pd <- pd[pd$Annotationlv2_step2 %in% cl,]
cells <- row.names(pd)
## 2. load ArchR project to obtain peak matrix-------------
proj <- loadArchRProject(path = "comATACx/")
proj <- proj[cells,]

#
print("Adding UMAP:")

proj <- addUMAP( proj,
                name = "UMAP_sub",
                reducedDims = "IterativeLSI_int",
                minDist = 0.2,
                metric = "euclidean",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
del <- getEmbedding(proj,"UMAP_sub")
colnames(del) <- c("iUMAP1","iUMAP2")
pd <- cbind(pd,del[row.names(pd),])

save(pd,file="~/SCENICOut/PN/pdPNatac.RData")
## subse tto make CistopicObj
#subset further, looks like 3k cells should be good
cells <- c()
cls <- unique(pd$Annotationlv2_step2)
for(x in cls){
  cells_sel <- row.names(pd)[pd$Annotationlv2_step2==x]
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

save(pd,file="~/SCENICOut/PN/pdPNatacdwn.RData")
cells <- row.names(pd)

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
  
projDir = "/home/SCENICOut/PN/Cistopic"
if not os.path.exists(projDir):
  os.makedirs(projDir)

os.chdir(projDir)
#
count_mat = sparse.csr_matrix(r.PM, dtype="int")
cn = r.a
peaks=r.peaks
#Create cisTopic object
path_to_blacklist='/home/extrafiles/hg38-blacklist.v2_noncanchr.bed'
print(" Creating pycistopic object from peak count matrix  ")

##path to fragment file not used
cistopic_obj = create_cistopic_object(fragment_matrix=count_mat,
                                      cell_names=cn, region_names= peaks,
                                      path_to_blacklist=path_to_blacklist)


print(" Adding metadata to cistopic object  ")

cistopic_obj.add_cell_data(r.pd)
print(cistopic_obj)
#Save
print(" Saving ...  ")
pickle.dump(cistopic_obj,open(os.path.join(projDir,'PN_cistopic_obj.pkl'), 'wb'))

exit
print("Finished creating cisTopic object for: Pontine ")
q()
