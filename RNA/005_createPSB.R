## this script creates pseudobulked clusters for the entire HB
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(rlist)
})
setwd("~/DATA/Diemana/scepro/")

load("combinedcountsHB.RData")

load("plotdataHBv2.RData")


com.sce <- com.sce[,row.names(plot.data)]

psb <- aggregateAcrossCells(com.sce, id=plot.data$NND_cl)
psb <- logNormCounts(Cluster)
save(Cluster, file ="HBNNDCl_psbANN.RData")

q()
### create cluster for a subset of cells
args <- commandArgs(trailing=TRUE)
load("combinedcountsHB.RData")

load("plotdataHBv2.RData")
##Here Annotation column could be level1 clustering or class 
plot.data <- plot.data[plot.data$Annotation==args,]

com.sce <- com.sce[,row.names(plot.data)]

psb <- aggregateAcrossCells(com.sce, id=plot.data$Annotation)
psb <- logNormCounts(Cluster)
save(Cluster, file =paste0("HB",args,"_psb.RData"))
q()