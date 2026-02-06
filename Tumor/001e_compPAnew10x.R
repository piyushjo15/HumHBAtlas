suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(dplyr)
  library(tidyr)
  library(rlist)
  
})
#load  -----
DIR_OUT="~/MBsnANA/HGGana/LGG_data/DJ_PA"
setwd(DIR_OUT)

# ## 1. loading and combining data --------
load("PA74scepro.RData")
load("PA87scepro.RData")
load("PA102scepro.RData")
load("PA158scepro.RData")

load("PA74dec.RData")
load("PA87dec.RData")
load("PA102dec.RData")
load("PA158dec.RData")

# ## combining as a list
all.sce <- list(PA74=PA74,PA87=PA87,PA102=PA102,PA158=PA158)
all.dec <- list(PA74=PA74.dec,
                PA87=PA87.dec,PA102=PA102.dec,
                PA158=PA158.dec)

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
save(plot.data, file="plotdataPAnew10x261124.RData") ## fixed ind sample clustering
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

load("HGNCsexchr.RData")
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
save(com.sce_cnnc, com.hvg, file="cosnormPAnew10x.RData")
rm(com.sce_cnnc)

## normalized across batches
norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce_rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce_rnnc) <- "logcounts"
mdt <- logcounts(com.sce_rnnc)
save(mdt, file ="lgcountsPAnew10x.RData")

##integrated plotdata
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsPAnew10x.RData")
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
