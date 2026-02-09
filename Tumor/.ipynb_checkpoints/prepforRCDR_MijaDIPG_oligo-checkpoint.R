## This script prepares input gene correlation matrix and gene list for graph 
## generation step for RECORDR analysis.The script generates this per sample
## i use pearson residual to obtain correlation
## Only focusing on Oligo component
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/p541i/MBsnANA/conda_env/pyforR/bin/python")
use_python("/home/p541i/MBsnANA/conda_env/pyforR/bin/python", required = TRUE)
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
})
options(future.globals.maxSize = 50 * 1024^3) 
##Directories 
DIR_IN="~/MBsnANA/HGGana/Mija_DIPG/"
setwd(DIR_IN)

DIR_OUT <- "~/MBsnANA/Diemana/scepro/RECORDR/"
sam <- commandArgs(trailingOnly = TRUE)
DIR_OUT <- paste0(DIR_OUT,sam)

print(paste0("Preparing generating input data for RECORDR for sample ",sam,"..."))

### 1.1 load expression and metadata ------
load("lgcountsMijaDIPG10x.RData")
load("plotdataMijaDIPG10x_ANN.RData")

## Removing ribosomal, mitochondiral and genes located on sex chromosomes
load("~/MBsnANA/Diemana/scepro/genesnoMTRibo.RData")
load("~/MBsnANA/Diemana/scepro/HGNCsexchr.RData")
genes <- genes[!(genes %in% SX)]

table(row.names(plot.data)==colnames(mdt))
## for v2
plot.data <- plot.data[plot.data$Batch==sam ,] #V2

cells <- readLines(paste0(DIR_OUT,"/",sam,"_cellv2_oligo_pr.txt"))
plot.data <- plot.data[cells,]
                   
## subset count matrix
mdt <- mdt[genes,row.names(plot.data)]
mdt[1:4,1:4]

## remove lowly expressed genes
genes <- row.names(mdt)
rsum <- rowSums(mdt>0)
rsum <-(rsum/dim(mdt)[2])*100
genes_sel <-genes[rsum>5]


### subsetting to only protein coding genes
pc <- readLines("~/DATA/HB/ann/gencdH38p13r37_PC.txt")
genes_sel <- genes_sel[genes_sel %in% pc]

mdt <- mdt[genes_sel,]
setwd(DIR_OUT)
repl_python()

import numpy as np
from scipy.sparse import csr_matrix, save_npz

mat = csr_matrix(r.mdt)
genes_sel=r.genes_sel
SAM=r.sam

save_npz(f"genemat_"+SAM+"_oligo_pc.npz", mat)
np.save(f"genes_oligo_pc.npy", genes_sel)
exit

