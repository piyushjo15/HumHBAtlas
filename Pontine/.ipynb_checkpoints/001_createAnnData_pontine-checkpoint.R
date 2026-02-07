## This scripts converts the RNA RData into python anData, 
library(reticulate)


Sys.setenv(RETICULATE_PYTHON = "/path/to/conda/envs/scenicplus/bin/python3.8")
use_python("/path/to/conda/envs/scenicplus/bin/python3.8")
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(scater)
})

DIR <- "~/MBsnANA/HBana/ATACana/Outs/"
DIR_RNA <- "~/MBsnANA/Diemana/scepro/"
DIR_SCENIC <- "~/MBsnANA/HBana/SCENICana/SCENICOut/PN/"

setwd(DIR)
## --------
##loading log normalized RNA data, no imputed cells in this sample
load(paste0(DIR_RNA,"pdPNdwn.RData"))
pd$Annotationlv2_step2 <- pd$Annotationlv2
load(paste0(DIR_SCENIC,"PN_HVG_pyscenic.RData"))
hvg <- hvg[1:2000]
hvg <- unique(c(hvg,heg,tfs_sel))

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

projDir = "/home/SCENICOut/PN/"

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
adata.write(os.path.join(projDir, 'adata_RNA.h5ad'), compression='gzip')

exit
q()
