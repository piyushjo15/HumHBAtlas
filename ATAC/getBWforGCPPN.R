
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
DIR <- "~/MBsnANA/HBana/ATACana/Outs/"

setwd(DIR)
#

proj <- loadArchRProject(path = "comATACx/")
load("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/plotdataHBATAC_210725.RData")

## oligo
cl_lv <- c("Neuroblasts_A:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_UBCP",
           "GCUBC:GCP_EM","GCUBC:GCP_PN")
keep <-  plot.data$Annotationlv2_step2 %in% cl_lv
plot.data<- plot.data[keep,]


cells <- row.names(plot.data)
proj <- proj[cells,]
proj$Annotationlv2_step2 <- plot.data$Annotationlv2_step2
getGroupBW(
  ArchRProj = proj,
  groupBy = "Annotationlv2_step2",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 2000,
  ceiling = 8,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
q()
