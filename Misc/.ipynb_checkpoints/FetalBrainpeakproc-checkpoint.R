## Processing peak data from Fetal brain study,
## identifying marker peaks for cluster annotaitons
library(tidyverse)

setwd("~/DeepHB/Outs/FetalBrain/")

### loading cluster metadata----
df <- read.delim("Cluster_metadata.txt")
df[1:4,1:4]
df$X <- NULL
df$Clusters <- paste0("Cl",df$Clusters)
##load marker peak matrix ------
peaks <- read.delim("All_peak_FetalBrain.500bp.strtmin1.bed",header = FALSE)
peaks$V2 <- peaks$V2+1
allpeaks <- row.names(peaks) <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
peakmat <- read.delim("Marker_peak_meta.txt")
peakmat$X <- NULL
colnames(peakmat) <-df$ClusterName
row.names(peakmat) <- allpeaks
peakmat[1:4,1:4]

cls <- df$ClusterName
markerpeaks <- list()
for(x in cls){
  keep <- peakmat[,x]>0
  markerpeaks[[x]] <- allpeaks[keep]
  rm(keep)
}
rm(x)
save(markerpeaks,peaks,df, file = "Allmarkerpeaks500bp.RData")
## generate bed file for marker peaks
load("Allmarkerpeaks500bp.RData")

##extracting all peaks
cls <- unique(df$ClusterName)
for(x in cls){
  sel_peaks <- markerpeaks[[x]]
  del <- peaks[sel_peaks,]
  del$V4 <- row.names(del)
  write.table(del, file = paste0("markerpeaks/",x,"_peaks.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              sep = "\t")
  rm(sel_peaks,del)
  
}

sel_class <- c("Oligo","Immune")
sel1 <- df$ClusterName[df$Class %in% sel_class]
sel2 <- grep("CB_GNP",df$ClusterName,value = TRUE)
sel3 <- grep("Purk",df$ClusterName,value = TRUE)
sel <- unique(c(sel1,sel2,sel3))

for(x in sel){
  sel_peaks <- markerpeaks[[x]]
  del <- peaks[sel_peaks,]
  del$V4 <- row.names(del)
  write.table(del, file = paste0("markerpeaks/",x,"_peaks.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              sep = "\t")
  rm(sel_peaks,del)
  
}








