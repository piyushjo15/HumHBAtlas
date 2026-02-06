
library(batchelor)
library(scran)
library(scater)
load("combinedcountsHB.RData")
assayNames(com.sce) <- "counts"

load("plotdataHBv2.RData")
keep <- row.names(plot.data.com)  %in% colnames(com.sce)
table(keep)
cells <- row.names(plot.data.com)[keep]
plot.data.com <- plot.data.com[cells,]
com.sce <- com.sce[,cells]
psb <- aggregateAcrossCells(com.sce, id=plot.data.com$Cluster)
psb <- logNormCounts(psb)
save(psb, file ="HBpsbANN.RData")
q()
save(del, file = "del.RData")
#failed due to zero counts for a gene
del2 <- round(counts(del))
del3 <- rowSums(del2)
x <- which(del3==min(del3))

del2 <- del2[-c(x),]
del3 <- SingleCellExperiment(list(counts=del2))
del3 <- logNormCounts(del3)
save(del3, file ="HBpsbANN.RData")
q()