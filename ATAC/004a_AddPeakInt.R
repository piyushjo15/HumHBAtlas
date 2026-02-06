## Calling peaks for the integrated ATAC data
suppressPackageStartupMessages({
  library(ArchR)
  library(scran)
  library(scater)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(tidyverse)
})

set.seed(456)
addArchRThreads(threads = 6)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/ATACana/Outs/"
setwd(DIR)
## 1. load data ----
proj <- loadArchRProject(path = "comATACx/") 

## loading preadjusted metadata
load("extrafiles/plotdataHB_ATACv6.RData")
ells <- row.names(plot.data.ATAC)

## subset cells
proj <- proj[cells,]

### Group coverages --------

#Creating PseudoBulk Replicates
proj <- addGroupCoverages(
  ArchRProj = proj,
  sampleRatio = 0.8,
  returnGroups = F,
  force = T,
  groupBy = "Clusters_int" , ## These are lv1 ATAC clusters, not label trasnfered clusteres
  minCells = 100,
  maxCells = 5000,
  minReplicates = 2,
  maxReplicates = 8,
  maxFragments = 100 * 10^6,
  useLabels= TRUE ## clusters are already defined per sample
)
print("Saving project after group coverage!")
saveArchRProject(proj, outputDirectory = "comATACx/")
#q()
# ###### Peak Calling ##########
#proj <- loadArchRProject(path = "comATACx/")

#Define path to MACS2
pathToMacs2 <- findMacs2()

#Peak Calling
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters_int",
  maxPeaks = 300000,##made it 300000
  pathToMacs2 = pathToMacs2,
  reproducibility = "2", ## ## 4 give errors
  excludeChr = c("chrY", "chrMT"),
  method = "q",
  cutOff = 0.01,
  extendSummits = 250,
  force = T
)
print("Saving project after adding peak matrix!")

saveArchRProject(proj, outputDirectory = "comATACx/")

PeakSet <- getPeakSet(proj)
save(PeakSet, file = "PeakCom/PeakSet.RData")
#q()
##------
#proj <- loadArchRProject(path = "comATACx/")
#load("PeakCom/PeakSet.RData")
keep <- PeakSet$Reproducibility>=4
PeakSet <- PeakSet[keep,]
save(PeakSet, file = "PeakCom/PeakSet_fil.RData")
proj <- addPeakSet(proj,peakSet = PeakSet,
                   force=TRUE)

#Adding Peak matrix
proj <- addPeakMatrix(proj,
                          ceiling = 10, ##changed it to 10
                          binarize = F,
                          force = T)
print("Saving project after adding peak matrix!")

saveArchRProject(proj, outputDirectory = "comATACx/")
q()

