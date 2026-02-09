##This script combines snRNA_seq data from PA samples
## The were intially processed similar to normal libraries
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(dplyr)
  library(tidyr)
  library(rlist)
  
})
#load  -----
DIR_OUT="~/PA"
setwd(DIR_OUT)

# ## 1. loading and combining data --------
load("PA01scepro.RData")
load("PA802scepro.RData")
load("PA03scepro.RData")

load("PA001dec.RData")
load("PA02dec.RData")
load("PA03dec.RData")

# ## combining as a list
all.sce <- list(PA001=PA001,PA02=PA02,PA03=PA03)
all.dec <- list(PA001=PA001.dec,PA02=PA02.dec,PA03=PA03.dec)

rm(list=ls(pattern="^PA"))
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
save(plot.data, file="plotdataPA.RData") ## fixed ind sample clustering
q()
load("comdatadec.RData")
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

save(com.hvg, file="topHVGPAnew10xnoRPMTSX.RData")
##cosine norm on HVG
com.sce_cnnc <- c()
for(x in sams){
  del <- logcounts(all.sce[[x]])
  del <- cosineNorm(del[com.hvg,])
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del)
}
save(com.sce_cnnc, com.hvg, file="cosnormPA.RData")
rm(com.sce_cnnc)

## normalized across batches
norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce_rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce_rnnc) <- "logcounts"
mdt <- logcounts(com.sce_rnnc)
save(mdt, file ="lgcountsPA.RData")

##integrated plotdata
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsPA.RData")
q()
## 3. HVG from tumor cells-----
sam <- read.delim("PAsams.txt",row.names = 1)
sam <- row.names(sam)

