## This script factorizes individual class+split using LDA into multiple ranks
## for meta-gene program analysis

library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")


suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(Matrix)
  library(batchelor)
})
setwd("~/MBsnANA/Diemana/scepro/")
args <- commandArgs(trailingOnly = TRUE)
clss <- args[1]
rnk <- args[2]
rnk <- as.integer(as.numeric(rnk))
Split <- "A"# or B
print(paste0("Identifying HVG for class: ",clss," and split: ",Split," and rank: ",rnk))
## 1. Processing raw counts for NMF run -------
#load raw counts
load("combinedcountsHB.RData")

#load plotdata and subset to selection
load("NMF_MP/plotdataHBforNMFMP.RData")
plot.data <- plot.data[plot.data$Class==clss & plot.data$SplitLDA==Split,]

#load HVG for class as calculated for scenic analysis
load(paste0("HVGs/",clss,"/",clss,"_HVGnTF.RData"))
if(length(hvg)>2000){
  hvg <- hvg[1:2000]
}
# downsampling by cluster
tb <- data.frame(table(plot.data$Cluster))
cl <- as.character(tb$Var1)[tb$Freq>3000]
if(length(cl)>0){
  keep <- plot.data$Annotationlv2 %in% cl
  table(keep)
  pds1 <- plot.data[keep, ]
  pds2 <- plot.data[!keep, ]
  pds3 <- c()
  rm(keep,tb)
  ##downsampling large annotations to max 3k cells
  for(y in cl){
    del <- pds1[pds1$Annotationlv2==y,]
    rn <- row.names(del)
    del2 <- del[sample(rn,3000),]
    pds3 <- rbind(pds3,del2)
    rm(del,rn,del2)
  }
  rm(y)
  plot.data <- rbind(pds2,pds3)
  rm(pds1,pds2,pds3)
}

# subset raw coutns
mdt <- t(round(counts(com.sce[hvg,row.names(plot.data)]))) 

## 2. python code for NMF ------
repl_python()
import sklearn.decomposition as sk
import numpy
import random

rnk=r.rnk
model = sk.LatentDirichletAllocation(n_components=rnk,random_state=0)
random.seed(456)
mdt=r.mdt.astype(int)
w = model.fit_transform(mdt)
h = model.components_

exit

rd_lda= py$w
H = py$h

row.names(rd_lda) <- row.names(plot.data)
colnames(rd_lda) <- row.names(H) <- paste0(clss,"_",Split,"_",rnk,"_",seq(dim(rd_lda)[2]))
colnames(H) <- hvg
save(rd_lda, plot.data, H, hvg, file = paste0("LDA_MP/",clss,"_",Split,"_",rnk,"_LDA.RData")) 
print("Finished!!")
q()
