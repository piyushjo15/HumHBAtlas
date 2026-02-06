## This script prepares input gene expression matrix and gene list for graph 
## generation step for RECORDR analysis.The script generates this per sample
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/p541i/MBsnANA/conda_env/pyforR/bin/python")
use_python("/home/p541i/MBsnANA/conda_env/pyforR/bin/python", required = TRUE)
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  
})

##Directories 
DIR_IN="~/MBsnANA/HGGana/LGG_data/DJ_PA/"
setwd(DIR_IN)

DIR_OUT <- "~/MBsnANA/Diemana/scepro/RECORDR/"
sam <- commandArgs(trailingOnly = TRUE)
DIR_OUT <- paste0(DIR_OUT,sam)

print(paste0("Preparing generating input data for RECORDR for sample ",sam,"..."))

### 1.1 load expression and metadata ------
load("lgcountsPAnew10x.RData")
##load metadata, combined
load("plotdataPAnew10x_ANN.RData")
table(row.names(plot.data)==colnames(mdt))

## Subset for a sample
#plot.data <- plot.data[plot.data$Batch==sam & plot.data$ANN1=="Tumor",] #v1

## Here subsetting for a sample and then extracting tumor cells annotated as
## OPC, Astro or Microglia
## for v2
plot.data <- plot.data[plot.data$Batch==sam ,] #V2
##select Astro, OPC and mciroglia
k1 <- grep("Astro",plot.data$ANN2)
k2 <- grep("OPC",plot.data$ANN2)
k3 <- grep("Microglia",plot.data$ANN3)
plot.data <- plot.data[c(k1,k2,k3),]

## This is just some extra step to remove cell annotation if the have less than
## 100 cells
## remove ANN which are less than 100 cells
tb <- data.frame(table(plot.data$ANN2))
tb <- as.character(tb$Var1)[tb$Freq>99]
plot.data <- plot.data[plot.data$ANN2%in%tb,]


## Within a sample sub-setting a tumor compartment to max 3k cells
## Most tumors have around 3-5k cells per compartment
cls <- unique(plot.data$ANN2)
ncells <- 3000
cells <- c()
set.seed(123)
for( x in cls){
  del_cells <- row.names(plot.data)[plot.data$ANN2==x]
  if(length(del_cells)>ncells){
    del_cells <- sample(del_cells,ncells)
  }
  cells <- c(cells,del_cells)
  rm(del_cells)
}
plot.data <- plot.data[cells,]

if(!dir.exists(DIR_OUT)){
  dir.create(DIR_OUT)
  
}
setwd(DIR_OUT)
write(cells, file = paste0(sam,"_cellv2.txt"))

## Here removing Ribosomal, mitochondrial and genes located in sex chromosomes
load("~/MBsnANA/Diemana/scepro/genesnoMTRibo.RData")
load("~/MBsnANA/Diemana/scepro/HGNCsexchr.RData")
genes <- genes[!(genes %in% SX)]

mdt <- mdt[genes,row.names(plot.data)]

## remove genes not expresse din at least 5% of cells
mdtx <- as(matrix(0,dim(mdt)[1],dim(mdt)[2]),"dgCMatrix")
mdtx[mdt>0] <-1
row.names(mdtx) <- row.names(mdt)
colnames(mdtx) <- colnames(mdt)
mdtx[1:4,1:4]

## Gene filtering step, a gene needs to be expressed in at least 5% of a cells
## in a tumor compartment
genes_sel <- c()
genes <- row.names(mdt)
for(x in cls){
  del_cells <- row.names(plot.data)[plot.data$ANN2==x]
  matx_del <-mdtx[,del_cells]
  rsum <- rowSums(matx_del)
  rsum <-(rsum/dim(matx_del)[2])*100
  genes_sel <- c(genes_sel,genes[rsum>5])
  rm(del_cells,matx_del,rsum)
}
genes_sel <- unique(genes_sel)

### subsetting to only protein coding genes
pc <- readLines("~/DATA/HB/ann/gencdH38p13r37_PC.txt")
genes_sel <- genes_sel[genes_sel %in% pc]

mdt <- mdt[genes_sel,]
rm(mdtx)


## find correlation matrix

repl_python()

import numpy as np
from scipy.sparse import csr_matrix, save_npz

mat = csr_matrix(r.mdt)
genes_sel=r.genes_sel
SAM=r.sam

save_npz("genemat_"+SAM+"_pcv2.npz", mat)
np.save("genes_pcv2.npy", genes_sel)
exit

##