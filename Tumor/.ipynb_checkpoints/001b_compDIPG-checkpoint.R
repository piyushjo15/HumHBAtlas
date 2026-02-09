##This script combines snRNA_seq data from DMG samples
## The were intially processed similar to normal libraries

suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(tidyverse)
  library(rlist)
  
})
#load  -----
DIR_OUT="~/DIPG/"
setwd(DIR_OUT)

# ## 1. loading and combining data --------
load("DMG01scepro.RData")
load("DMG02cepro.RData")
load("DMG03scepro.RData")
load("DMG04cepro.RData")
load("DMG05scepro.RData")

load("DMG01dec.RData")
load("I007_107dec.RData")
load("DMG03dec.RData")
load("I027_024dec.RData")
load("DMG05dec.RData")

# ## combining as a list
all.sce <- list(DMG01=DMG01,I007_107=I007_107,DMG03=DMG03,
                I027_024=I027_024,DMG05=DMG05)
all.dec <- list(DMG01=DMG01.dec,
                I007_107=I007_107.dec,DMG03=DMG03.dec,
                I027_024=I027_024.dec,DMG05=DMG05.dec)

rm(list=ls(pattern="^I0"))
universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce <- lapply(all.sce, "[", i=universe,)
save(all.sce, file="comdata.RData")

universe2 <- Reduce(intersect, lapply(all.dec, rownames))
all.dec <- lapply(all.dec, "[", i=universe2,)
save(all.dec, file="comdatadec.RData")

## 2. combining for getting an integrated gene exp mat ----
#load("comdata.RData")

#plot data
l <- length(all.sce)

plot.data <- c()
for(i in 1:l){
  sce <- all.sce[[i]]
  del <- data.frame(colData(sce))
  del2 <- data.frame(reducedDim(sce, "UMAP"))
  colnames(del2) <- c("iUMAP1", "iUMAP2")
  del <- cbind(del,del2)
  plot.data <- rbind(plot.data,del)
  rm(del,del2)
}
save(plot.data, file="plotdataDIPG.RData") 
#q()
#load("comdatadec.RData")
## HVG
sams <- names(all.dec)
com.hvg <- com.hvg2 <- c()
for(x in sams){
  del <- getTopHVGs(all.dec[[x]])
  com.hvg <- c(com.hvg,del[1:1500])
  com.hvg2 <- c(com.hvg2,del[1:3000])
}
rm(x)
com.hvg <- com.hvg[!duplicated(com.hvg)]

rm(all.dec)

load("~/extrafiles/HGNCsexchr.RData")
com.hvg <- com.hvg[!(com.hvg %in% SX)]

save(com.hvg, file="topHVGDIPGnoRPMTSX.RData")
##cosine norm on HVG
com.sce_cnnc <- c()
for(x in sams){
  del <- logcounts(all.sce[[x]])
  del <- cosineNorm(del[com.hvg,])
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del)
}
save(com.sce_cnnc, com.hvg, file="cosnormIPG.RData")
rm(com.sce_cnnc)

## normalized across batches
norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce_rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce_rnnc) <- "logcounts"
mdt <- logcounts(com.sce_rnnc)
save(mdt, file ="lgcountsDIPG.RData")

##integrated plotdata
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsDIPG.RData")
q()

