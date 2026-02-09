## This script prepares input gene correlation matrix and gene list for graph 
## generation step for RECORDR analysis.The script generates this per sample
## i use pearson residual to obtain correlation

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(glmGamPoi)
  library(Seurat)
})
options(future.globals.maxSize = 50 * 1024^3) 
##Directories 
DIR_IN="~/PA/"
setwd(DIR_IN)

DIR_OUT <- "~/RECORDR/"
sam <- commandArgs(trailingOnly = TRUE)
DIR_OUT <- paste0(DIR_OUT,sam)

print(paste0("Preparing generating input data for RECORDR for sample ",sam,"..."))

### 1.1 load expression and metadata ------
load("combinedcountsPA.RData")
mdt <- counts(com.sce)
rm(com.sce)
##load metadata, combined
load("plotdataPA_ANN.RData")
table(row.names(plot.data)==colnames(mdt))

## Removing ribosomal, mitochondiral and genes located on sex chromosomes
load("~/extrafiles/genesnoMTRibo.RData")
load("~/extrafiles/HGNCsexchr.RData")
genes <- genes[!(genes %in% SX)]

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
## subset count matrix
mdt <- mdt[genes,row.names(plot.data)]
mdt[1:4,1:4]

## remove lowly expressed genes
genes_sel <- c()
genes <- row.names(mdt)
for(x in cls){
  del_cells <- row.names(plot.data)[plot.data$ANN2==x]
  mdtx_del <- mdt[,del_cells]
  rsum <- rowSums(mdtx_del>0)
  rsum <-(rsum/dim(mdtx_del)[2])*100
  genes_sel <- c(genes_sel,genes[rsum>5])
  rm(del_cells,mdtx_del,rsum)
}
genes_sel <- unique(genes_sel)

### subsetting to only protein coding genes
pc <- readLines("~/extrafiles/gencdH38p13r37_PC.txt")
genes_sel <- genes_sel[genes_sel %in% pc]

mdt <- mdt[genes_sel,]
## using Seurat SCTtransform to find pearson residuals for correlation
sct <- CreateSeuratObject(counts = mdt,meta.data = plot.data)
sct <- SCTransform(sct,ncells = 3000, 
                   vst.flavor = "v2",return.only.var.genes = FALSE)
R <- sct@assays$SCT$scale.data
dim(R)
R[1:3,1:3]
R <- R[genes_sel,row.names(plot.data)]
rm(mdt,sct)
## find correlation matrix

repl_python()

import numpy as np
from scipy.sparse import csr_matrix, save_npz

mat = csr_matrix(r.R)
genes_sel=r.genes_sel
SAM=r.sam

save_npz(f"genemat_"+SAM+"_pc.npz", mat)
np.save(f"genes_pc.npy", genes_sel)
exit

