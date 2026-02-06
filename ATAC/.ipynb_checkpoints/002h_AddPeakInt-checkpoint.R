## getting peaks called on cluster_subset
suppressPackageStartupMessages({
  library(ArchR)
  library(scran)
  library(scater)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(dplyr)
  library(tidyr)
})

set.seed(456)
addArchRThreads(threads = 6)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/MBsnANA/HBana/ATACana/Outs/"
DIRS <- "~/MBsnANA/HBana/ATACana/HBATAC_Scripts/"

setwd(DIR)
#
# ## 1. load data ----
# proj <- loadArchRProject(path = "comATACx/") ## this was created by me before by copying comATAC
# 
# ## loading preadjusted metadata
# load("extrafiles/plotdataHB_ATACv5.RData")
# ## remove subclusters
# torem <- readLines("extrafiles/Clusters_subset_toremove.txt")
# keep <- plot.data.ATAC$Clusters_subset %in% torem
# table(keep)
# cells <- row.names(plot.data.ATAC)[!keep]
# 
# ## subset cells
# proj <- proj[cells,]
# 
# ### Group coverages --------
# 
# #Creating PseudoBulk Replicates
# proj <- addGroupCoverages(
#   ArchRProj = proj,
#   sampleRatio = 0.8,
#   returnGroups = F,
#   force = T,
#   groupBy = "Clusters_int" ,
#   minCells = 100,
#   maxCells = 5000,
#   minReplicates = 2,
#   maxReplicates = 8,
#   maxFragments = 100 * 10^6,
#   useLabels= TRUE ## clusters are already defined per sample
# )
# print("Saving project after group coverage!")
# saveArchRProject(proj, outputDirectory = "comATACx/")
# q()
# # ###### Peak Calling ##########
# proj <- loadArchRProject(path = "comATACx/")
# 
# #Define path to MACS2
# pathToMacs2 <- findMacs2()
# 
# #Peak Calling
# proj <- addReproduciblePeakSet(
#   ArchRProj = proj,
#   groupBy = "Clusters_int",
#   maxPeaks = 300000,##made it 300000
#   pathToMacs2 = pathToMacs2,
#   reproducibility = "2", ## ## 4 give errors
#   excludeChr = c("chrY", "chrMT"),
#   method = "q",
#   cutOff = 0.01,
#   extendSummits = 250,
#   force = T
# )
# print("Saving project after adding peak matrix!")
# 
# saveArchRProject(proj, outputDirectory = "comATACx/")
# 
# PeakSet <- getPeakSet(proj)
# save(PeakSet, file = "comATACx/PeakSet.RData")
# ##------
# proj <- loadArchRProject(path = "comATACx/")
# load("comATACx/PeakSet.RData")
# keep <- PeakSet$Reproducibility>=4
# PeakSet <- PeakSet[keep,]
# save(PeakSet, file = "comATACx/PeakSet_fil.RData")
# proj <- addPeakSet(proj,peakSet = PeakSet,
#                    force=TRUE)
# 
# #Adding Peak matrix
# proj <- addPeakMatrix(proj,
#                           ceiling = 10, ##changed it to 10
#                           binarize = F,
#                           force = T)
# print("Saving project after adding peak matrix!")
# 
# saveArchRProject(proj, outputDirectory = "comATACx/")
# q()
# ## Converting PeakSet to bed -------
# load("comATACx/PeakSet_fil.RData")
# 
# PS <- data.frame(PeakSet) 
# row.names(PS) <- paste0(PS$seqnames,":",PS$start,"-",PS$end)
# head(PS)
# bed <- PS[,c(1:3)]
# bed$X <- row.names(PS)
# write.table(bed, file="comATACx/AllPeaks_robust.bed", col.names = FALSE, row.names = FALSE,quote = FALSE,sep="\t")
# ##then i remvoed HG38 blacklisted region and noncanonical chromosomes from it
# q()
# ##adding cluster info
# cl <- PS$GroupReplicate
# cl <-gsub("._."," ",cl,fixed = TRUE)
# cl <- data.frame(X=cl)
# 
# cl2 <- cl %>% separate(X,c("A","B"), sep=" ")
# head(cl2)
# PS$Cluster <- cl2$A
# head(PS)
# rm(cl,cl2)
# table(row.names(PS)==row.names(region_names))
# ##get peak matrix-----
# proj <- loadArchRProject(path = "comATACx/")
# 
# PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
# save(PeakMatrix, file = "comATACx/PeakMatrix.RData")
# q()
## get marker peaks -----
proj <- loadArchRProject(path = "comATACx/")
load("extrafiles/plotdataHB_ATACv6_290125.RData")
cells1 <- getCellNames(proj)
cells2 <- row.names(plot.data.ATAC)
cells1 <- cells1[cells1 %in% cells2]
plot.data.ATAC <- plot.data.ATAC[cells1,]
proj <- proj[cells1,]
proj$Ann <-plot.data.ATAC$Class
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Ann",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 5000,k = 300
)
save(markersPeaks, file = "extrafiles/markerpeaks_Class.RData")
q()
load("extrafiles/markerpeaks_Class.RData")
allpeaks <- read.delim("PeakCom/AllPeaks_robust.fil.hg19v2.bed",header = FALSE)
rn <- row.names(allpeaks) <- allpeaks$V4
head(allpeaks)

markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.01 & (Log2FC >=2 | AUC >= 0.52)")
lengths(markerList)

sc <- names(markerList)
for(x in sc){
  del <- data.frame(markerList[[x]])
  del <- del[,c(1,3,4)]
  id <- paste0(del$seqnames,":",del$start,"-",del$end)
  del$peak <- id
  write.table(del,paste0("SNP_pGlioma/",x,"_CREs.bed"),
              sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
  # ##hg19
  # rnx <- rn[rn %in% id]
  # ## markerv4 is union of wilcox and T, hg19
  # write.table(allpeaks[rnx,],paste0("GWAS/markersv4/",x,"_CREs.bed"),
  #             sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

  rm(del)
}
rm(x)
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR < 0.01 & (Log2FC >=2 | AUC >= 0.52)",
  transpose = TRUE
)
mdt <- heatmapPeaks@matrix
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#add background peaks----
proj <- loadArchRProject(path = "comATACx/")
proj <- addBgdPeaks(proj, force = T)

load("comATACx/PeakSet.RData")
head(PeakSet)
table(PeakSet$Reproducibility)
