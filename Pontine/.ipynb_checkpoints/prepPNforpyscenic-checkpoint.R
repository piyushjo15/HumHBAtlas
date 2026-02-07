## Preparing a loom file for each Class
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(Matrix)
  library(tidyverse)
  
})

#Directories
DIR_RNA <-"~/RNA/"
DIR_OUT <-"~/SCENICOut/PN"
setwd(DIR_OUT)

## 1. --------
#loading counts file
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
#load hvg data
load("PN_HVG_pyscenic.RData")
hvg <- hvg[1:2000]
hvg <- unique(c(hvg,heg,tfs_sel))
#load plotdata
load(paste0(DIR_RNA,"pdPNdwn.RData"))

library(SCENIC)
library(SCopeLoomR)
##loom file
print("Creating data for loom file")
exprMat <- round(counts(com.sce)[hvg,row.names(pd)])
cellInfo <- pd
SCE <- build_loom("PN.loom",dgem=exprMat)
SCE <- add_cell_annotation(SCE, cellInfo)
close_loom(SCE)
print("Bye!")

q()
