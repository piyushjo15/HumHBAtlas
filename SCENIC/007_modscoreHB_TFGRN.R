## This script obtained module score for all the TF-GRNs identified across classes
## To run the analysis without overloading memory requirements, the scores are calculated
## for GRNS per class and then merged
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(tidyverse)
  library(AUCell)
  library(BiocParallel)
})
##this fucntion has been provided in other 'modscore' scripts
## removed here for space
addmodscore <- function()
setwd("~/RNA/")
Cls <- commandArgs(trailing=TRUE)
##1. load required data -------
### 1.1 load gene-sets ------
# ##TF-GRNs
load("~/SCENIC/AllGRNs.RData")
grns <- row.names(all_tfs) <- paste0(all_tfs$TF,"_",all_tfs$Class)
sel_grns <- grns[all_tfs$Class==Cls]
GSEA <- list()
for(x in sel_grns){
  GSEA[[x]] <- all_grns[[x]]
}
### 1.2 load expression and metadata ------

## load expression data
## this is compiled log transformed expression data for the entire hindbrain
load("lgcountsHB.RData")

##load metadata, combined
load("plotdataHBv2.RData")
plot.data <-plot.data[plot.data$Class=Cls,]
mdt <- mdt[,row.names(plot.data)]

## 2. module score----
scr <- addmodscore(mdt,features = GSEA)
colnames(scr) <- names(GSEA)
row.names(scr) <- row.names(plot.data)

save(scr, plot.data, file = paste0("AUC_MS/pdHB_",Cls,"_MS.RData"))

q()

## now merging all the GRNs
zscore_c <- function(x){
  sdx <- sd(x)
  mx <- mean(x)
  ot <- boxplot.stats(x)$out
  q <- boxplot.stats(x)$stats
  keep <- (x < q[1])
  x[keep] <- q[1]
  keep <- (x > q[5])
  x[keep] <- q[5]
  x2 <- x-mx
  z <- x2/sdx
  return(z)
}
Cls <- readLines("~/RNA/RNA_Class.txt")


scr_com <- c()
scr_com_z <- c()

for(x in Cls){
  load(paste0("AUC_MS/pdHB_",x,"_MS.RData"))
  scrx <- apply(scr, 2, zscore_c)
  scr_com <- cbind(scr_com,as.matrix(scr))
  scr_com_z <- cbind(scr_com_z,scrx)
  rm(scr)
}
rm(x)
save(scr_com,scr_com_z, plot.data, file = "AUC_MS/pdHB_AllGRNs_MSz.RData") 
q()