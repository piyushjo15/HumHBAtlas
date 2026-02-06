suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(tidyverse)
  library(rlist)
  
})
#load  -----
DIR_OUT="~/MBsnANA/HGGana/Mija_DIPG/"
setwd(DIR_OUT)

# ## 1. loading and combining data --------
load("I007_070scepro.RData")
load("I007_107scepro.RData")
load("I023_009scepro.RData")
load("I027_024scepro.RData")
load("I054_021scepro.RData")

load("I007_070dec.RData")
load("I007_107dec.RData")
load("I023_009dec.RData")
load("I027_024dec.RData")
load("I054_021dec.RData")

# ## combining as a list
all.sce <- list(I007_070=I007_070,I007_107=I007_107,I023_009=I023_009,
                I027_024=I027_024,I054_021=I054_021)
all.dec <- list(I007_070=I007_070.dec,
                I007_107=I007_107.dec,I023_009=I023_009.dec,
                I027_024=I027_024.dec,I054_021=I054_021.dec)

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
save(plot.data, file="plotdataMijaDIPG10x170625.RData") ## fixed ind sample clustering
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

load("~/MBsnANA/Diemana/scepro/HGNCsexchr.RData")
com.hvg <- com.hvg[!(com.hvg %in% SX)]

save(com.hvg, file="topHVGMijaDIPG10xnoRPMTSX.RData")
##cosine norm on HVG
com.sce_cnnc <- c()
for(x in sams){
  del <- logcounts(all.sce[[x]])
  del <- cosineNorm(del[com.hvg,])
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del)
}
save(com.sce_cnnc, com.hvg, file="cosnormMijaDIPG10x.RData")
rm(com.sce_cnnc)

## normalized across batches
norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce_rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce_rnnc) <- "logcounts"
mdt <- logcounts(com.sce_rnnc)
save(mdt, file ="lgcountsMijaDIPG10x.RData")

##integrated plotdata
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsMijaDIPG10x.RData")
q()
## 3. HVG from tumor cells-----
sam <- read.delim("PAnew10x.txt",row.names = 1)
sam <- row.names(sam)

##from combined plot data extract tumor cells (oligo+Astro)
load("plotdataPAnew10x_ANN.RData")
load("~/MBnewana/genesnoRPMT.RData")
load("HGNCsexchr.RData")

plot.data <- plot.data[plot.data$ANN1=="Tumor",]
com.hvg <- c()
sams <- unique(plot.data$Batch)
for(x in sams){
  sce <- get(load(paste0(x,"scepro.RData")))
  pd <- plot.data[plot.data$Batch==x,]
  rn <- row.names(pd)
  sce <- sce[genes,rn]
  sce.dec <- modelGeneVar(sce)
  top.hvgs <- getTopHVGs(sce.dec, n = 5000)
  top.hvgs <- top.hvgs[!(top.hvgs %in% SX )]
  com.hvg <- c(com.hvg,top.hvgs[1:2000])
  rm(sce,pd,rn,top.hvgs)
  rm(list=ls(pattern="^PA"))
}
rm(x)
com.hvg <- com.hvg[!duplicated(com.hvg)]
save(com.hvg, file="topHVGPAnew10xnoRPMTSX_TUM.RData")
