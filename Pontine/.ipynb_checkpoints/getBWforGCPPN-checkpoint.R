
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
})


set.seed(456)
addArchRThreads(threads = 8)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~ATACana/Outs/"

setwd(DIR)
#

proj <- loadArchRProject(path = "comATACx/")
load("plotdataHBATACv2.RData")

## oligo
cl_lv <- c("Neuroblasts_I:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_UBCP",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
keep <-  plot.data$Cluster_step2 %in% cl_lv
plot.data<- plot.data[keep,]


cells <- row.names(plot.data)
proj <- proj[cells,]
proj$Cluster_step2 <- plot.data$Cluster_step2
getGroupBW(
  ArchRProj = proj,
  groupBy = "Cluster_step2",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 2000,
  ceiling = 8,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
q()
