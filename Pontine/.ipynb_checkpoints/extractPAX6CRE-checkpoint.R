## Peak count matrix for HOX bdining sites on PAX6 locus
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ComplexHeatmap)
  library(scater)
  library(scran)
  library(AUCell)
  
})

DIR_ATAC <- "~/ATACana/Outs/"
DIR_GRN <- "~/SCENIC/"

df <- read.delim(paste0(DIR_GRN,"PN/metadata/eRegulons_fix.txt"))

hox <- unique(df$TF)
hox <- grep("HOX",hox,value = TRUE)
del <- df[df$TF %in% hox,]
del <- del[del$Gene=="PAX6",]
cres <- sort(unique(del$Region))
sel_cls <- c("Neuroblasts_I:A1_AES","Glioblasts:Progenitor_AES",
             "GCUBC:GCP_PrN","GCUBC:GCP_PN")

load("plotdataHBATACv2.RData")
plot.data <- plot.data[plot.data$Cluster_step2%in%sel_cls,]

cells <- row.names(plot.data)
## process subsetted peak matrix------------
proj <- loadArchRProject(path = paste0(DIR_ATAC,"comATACx/"))

##extract PeakMatrix for selected cells
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)

## get Matrix and subset to peaks
region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)

PM <- assay(PeakMatrix)
row.names(PM) <- row.names(region_names)
PM[1:4,1:4]
PM <- PM[,cells] ###very important to make sure cell order is same


### generate psuedobulk
mdt <- SingleCellExperiment(list(counts=PM))
psb <- aggregateAcrossCells(mdt,plot.data$Cluster_step2)
psb <- logNormCounts(psb)
mdt <- logcounts(psb) ## for Fig 3 i have used this
mdt[1:4,1:4]
mdt <- mdt[cres,]
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
mdt2 <- t(apply(mdt,1,range_01))
mdt2[1:4,1:4]
save(mdt2,mdt, file="PAX6_HOX_CRE.RData")
q()
mc <- colorRampPalette(c("white","grey75","black"))(101)
mb <- (seq(1:101)-1)/100