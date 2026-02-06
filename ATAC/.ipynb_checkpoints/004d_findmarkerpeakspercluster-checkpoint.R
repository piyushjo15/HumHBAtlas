## finding marker peaks per cluster by comparing to clusters within that class
suppressPackageStartupMessages({
  library(ArchR)

})
setwd("~/MBsnANA/HBana/ATACana/Outs/")
set.seed(456)
cls <- commandArgs(trailing=TRUE)
addArchRThreads(8)
print(paste0("finding marker CREs per cluster for class :",cls))
## part 1. load all the data ------------
### 1.1 load plotdata and subset to max 200 cells
load("plotdataHB_ATACv6.RData")
tb <- data.frame(table(plot.data.ATAC$Cluster_step2))
tb <- as.character(tb$Var1)[tb$Freq>49]
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$Cluster_step2 %in% tb,]
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$Class==cls,]
proj <- loadArchRProject(path = "comATACx/")

cells1 <- getCellNames(proj)
cells2 <- row.names(plot.data)
cells1 <- cells1[cells1 %in% cells2]
plot.data.ATAC <- plot.data.ATAC[cells1,]
proj <- proj[cells1,]
table(row.names(plot.data.ATAC)==getCellNames(proj))

proj$Ann <-plot.data.ATAC$Cluster_step2
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Ann",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 1000
)
save(markersPeaks, file = paste0("PeakCom/markerpeaks_cluster_",cls,".RData"))
q()

## find marker CRE per cluster
cls <- readLines("~/RNA/RNA_Class.txt")
all_peaks <- read.delim("PeakCom/AllPeaks_robust.fil.sorted.bed",header = FALSE,row.names = 4)
colnames(all_peaks) <- c("chr","start","end")
all_peaks$peak <- row.names(all_peaks)

marker_peaks <- c()
for(x in cls){
  load(paste0("ATAC/markerpeaks_cluster_",x,".RData"))
  markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.05 & (Log2FC >=2 | AUC >= 0.55)")
  clusters <- names(markerList)
  npeaks <- lengths(markerList)
  
  for(y in clusters){
    peaklen <- npeaks[y]
    if(peaklen>0){
      del <- data.frame(markerList[[y]])
      del <- del[,c(1,3,4)]
      colnames(del) <- c("chr","start","end")
      id <- paste0(del$chr,":",del$start,"-",del$end)
      del$peak <- id
      del$Cluster <- y
      if(peaklen>500){
        del <- del[1:500,]
      }
      marker_peaks <- rbind(marker_peaks,del)
      rm(del)
      
    }
    
  }
  rm(clusters,y,markerList,markersPeaks)
}
rm(x)


save(marker_peaks, file="PeakCom/markerpeaks_cluster_com.RData")

## add it to all peaks
library(tidyverse)

load("PeakCom/markerpeaks_cluster_com.RData")
marker_peaks <- marker_peaks[!is.na(marker_peaks$Cluster),]
peaks_sel <- unique(marker_peaks$peak)

res <- marker_peaks %>%
  filter(peak %in% peaks_sel) %>%
  group_by(peak) %>%
  summarise(Cluster = paste(Cluster, collapse = ","), .groups = "drop")

all_peaks <- left_join(all_peaks, res, by = c("peak"))
load("PeakCom/markerpeaks_class_com.RData")

peaks_sel <- unique(marker_peaks$peak)
res <- marker_peaks %>%
  filter(peak %in% peaks_sel) %>%
  group_by(peak) %>%
  summarise(Class = paste(Class, collapse = ","), .groups = "drop")

all_peaks <- left_join(all_peaks, res, by = c("peak"))
save(all_peaks, file="PeakCom/all_peaks_with_markers.RData")
write.table(all_peaks, file="PeakCom/AllPeaks_robust_ann.fil.sorted.bed",sep="\t",row.names=FALSE,quote=FALSE)

