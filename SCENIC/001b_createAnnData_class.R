## This scripts converts the RNA RData into python anData, that will be used as 
## SCENIC+ input. 
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/path/to/conda/envs/scenicplus/bin/python3.8")
use_python("/path/to/conda/envs/scenicplus/bin/python3.8")
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(scater)
})

DIR_ATAC <- "~/ATACana/Outs/"
DIR_RNA <- "~/RNA/"
DIR_OUT<- "~/SCENICOut/Class/"

setwd(DIR_ATAC)
Cls <- commandArgs(trailingOnly = TRUE)
Splt <- "A" ## "A" or "B"
print(paste0("Creating annData object for RNA data of class: ",Cls," and split: ",Splt))
## --------

#load plotdata and subset to selection
load(paste0(DIR_RNA,"plotdataHBv2.RData"))
pds <- plot.data[plot.data$Class==Cls & plot.data$SplitLDA==Splt,]
rm(plot.data)
# downsampling
tb <- data.frame(table(pds$Cluster))
cl <- as.character(tb$Var1)[tb$Freq>3000]
keep <- pds$Cluster %in% cl
table(keep)
pds1 <- pds[keep, ]
pds2 <- pds[!keep, ]
pds3 <- c()
rm(keep,tb)
##downsampling large annotations to max 3k cells
for(y in cl){
  del <- pds1[pds1$Cluster==y,]
  rn <- row.names(del)
  del2 <- del[sample(rn,3000),]
  pds3 <- rbind(pds3,del2)
  rm(del,rn,del2)
}
rm(y)
pd <- rbind(pds2,pds3)
rm(pds,pds1,pds2,pds3,cl)
## this is done as the 'Cluster' labels in ATAC data are stored in 
## Cluster_step2 column
pd$Cluster_step2 <-  pd$Cluster 
## for Neurons_B_B remove 
if(Cls=="Neurons_II"){
  clx <- c("Neurons_II:ISO","Neurons_II:LHX2_PRDM8")
  
  keep <- pd$Cluster %in% clx
  pd <- pd[!keep,]
}

#load HVG for class as calculated for scenic analysis
load(paste0("~/HVGs/",Cls,"/",Cls,"_HVGnTF.RData"))
if(length(hvg)>2000){
  hvg <- hvg[1:2000]
  
}
hvg <- unique(c(hvg,heg,tfs_sel))

##loading log normalized RNA data
load(paste0(DIR_RNA,"lgcountsHB.RData"))

mdt <- mdt[hvg,row.names(pd)]

#save the barcodes
barcodes<-data.frame(colnames(mdt))
colnames(barcodes)<-'Barcode'

#Save the gene names
genes<-data.frame(rownames(mdt))
colnames(genes)<-'Gene'


######### Save the data to python ############
repl_python()

import scanpy as sc
import pandas as pd
from scipy import io
import os
from scipy import sparse

DIR = '/path/to/SCENICOut/Class/'
Cls=r.Cls
Splt=r.Splt
projDir = os.path.join(DIR,Cls)

#log count Matrix
counts = sparse.csr_matrix(r.mdt)
lgcounts = r.mdt
#load metadata
barcodes=r.barcodes
genes=r.genes

adata=sc.AnnData(counts.T, dtype=counts.T.dtype)

adata.raw = adata
adata.layers["log1p"] = lgcounts.T
adata.obs_names=barcodes['Barcode'].values
adata.var_names=genes['Gene'].values

#import metadata
adata.obs=r.pd

#save adata object
adata.write(os.path.join(projDir,Cls+'_'+Splt+'_RNA.h5ad'), compression='gzip')

exit
print(paste0("Finished creating annData object for RNA data of class: ",Cls," and split: ",Splt))

suppressPackageStartupMessages({
  library(SCENIC)
  library(SCopeLoomR)
  library(scater)
})
rm(mdt)
## here I also create loom file for pyscenic run
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
exprMat <- round(counts(com.sce)[hvg,row.names(pd)])  
cellInfo <- pd

SCE <- build_loom(paste0(Cls,"/",Cls,"_",Splt,".loom"),dgem=exprMat)
SCE <- add_cell_annotation(SCE, cellInfo)
close_loom(SCE)
print(paste0("finished creating the loom file...",Cls,"/",Cls,"_",Splt,".loom"))


q()
