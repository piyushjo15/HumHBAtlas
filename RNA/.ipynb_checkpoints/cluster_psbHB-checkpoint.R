
library(batchelor)
library(scran)
library(scater)
#library(scuttle)
#args <- "Neurons"
setwd("~/DATA/Diemana/scepro/")

load("combinedcountsHB.RData")

load("plotdataHBv2.RData")


com.sce <- com.sce[,row.names(plot.data)]

psb <- aggregateAcrossCells(com.sce, id=plot.data$NND_cl)
psb <- logNormCounts(Cluster)
save(Cluster, file ="HBNNDCl_psbANN.RData")

q()
