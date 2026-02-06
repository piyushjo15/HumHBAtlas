## This script prepares input gene correlation matrix and gene list for graph 
## generation step for RECORDR analysis.The script generates this per sample
## i use pearson residual to obtain correlation

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/p541i/MBsnANA/conda_env/pyforR/bin/python")
use_python("/home/p541i/MBsnANA/conda_env/pyforR/bin/python", required = TRUE)
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(glmGamPoi)
  library(Seurat)
})
options(future.globals.maxSize = 50 * 1024^3) 
##Directories 
##Directories 
DIR_IN="~/MBsnANA/HGGana/PFA/"
setwd(DIR_IN)

DIR_OUT <- "~/MBsnANA/Diemana/scepro/RECORDR/"
sam <- commandArgs(trailingOnly = TRUE)
DIR_OUT <- paste0(DIR_OUT,sam)

print(paste0("Preparing generating input data for RECORDR for sample ",sam,"..."))

### 1.1 load expression and metadata ------
load("combinedcountsGJPFA.RData")
mdt <- counts(com.sce)
rm(com.sce)
##load metadata, combined
load("plotdataGJPFA_121124v2_ANN.RData")
table(row.names(plot.data)==colnames(mdt))

## for Gojo I have to find a subset of genes
load("genesnoRPMT_Gojo.RData")
genes1 <- genes
load("~/MBsnANA/Diemana/scepro/genesnoMTRibo.RData") ## For Gojo PFA I have to subset to tumor genes
genes2 <- genes
genes <- intersect(genes1, genes2)


plot.data <- plot.data[plot.data$Batch==sam ,] #V2
##select Ependymal, mciroglia
k1 <- grep("Ependymal",plot.data$ANN2)
k2 <- grep("Microglia",plot.data$ANN2)
plot.data <- plot.data[c(k1,k2),]

## remove ANN which are less than 100 cells
tb <- data.frame(table(plot.data$ANN2))
tb <- as.character(tb$Var1)[tb$Freq>99]
plot.data <- plot.data[plot.data$ANN2%in%tb,]
#### not subsettign anythng
if(!dir.exists(DIR_OUT)){
  dir.create(DIR_OUT)
  
}
setwd(DIR_OUT)

mdt <- mdt[genes,row.names(plot.data)]

## remove genes not expresse din at least 5% of cells in any annotion
genes_sel <- c()
genes <- row.names(mdt)
cls <- unique(plot.data$ANN2)
for(x in cls){
  del_cells <- row.names(plot.data)[plot.data$ANN2==x]
  matx_del <-mdt[,del_cells]
  rsum <- rowSums(matx_del>0)
  rsum <-(rsum/dim(matx_del)[2])*100
  genes_sel <- c(genes_sel,genes[rsum>5])
  rm(del_cells,matx_del,rsum)
}
genes_sel <- unique(genes_sel)
pc <- readLines("~/DATA/HB/ann/gencdH38p13r37_PC.txt")
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

save_npz(f"genemat_"+SAM+"_pc_pr.npz", mat)
np.save(f"genes_"+SAM+"_pc_pr.npy", genes_sel)
exit

