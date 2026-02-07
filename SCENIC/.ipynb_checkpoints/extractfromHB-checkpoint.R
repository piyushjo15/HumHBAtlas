##This script is to extract data from Hindbrain atlas
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
})

DIR_HB <-"~/DATA/Diemana/scepro/"
##set working directory for SCENIC data analysis
setwd("~/MBsnANA/HBana/SCENICana/SCENICOut/")
# load data and process ----
## Loading HB ref 
load(paste0(DIR_HB,"lgcountsHB.RData"))
load(paste0(DIR_HB,"plotdata_NNDCS_50k_180823fix.RData")) 
load(paste0(DIR_HB,"genesnocc.RData"))
#extracting Purkinje Matrue and GCP cells for analysis----
##using NND_cluster for now
#dvs <- c("NNDCS_145") ##selecting 
##using annotations
dvs <- c("X","Neuroblasts","Neuroblasts2","Neuroblasts3",
         "Neurons2","Neurons3","Purkinje") ##selecting 
keep <- plot.data$Annotation %in% dvs
table(keep)
plot.data <- plot.data[keep,]
mdt <- mdt[,row.names(plot.data)]
dec <- modelGeneVar(mdt)
common.genes <- getTopHVGs(dec)
common.genes <- common.genes[common.genes %in% genes]

save(mdt,plot.data, file = "HUMHB_NB_Neurons.RData")
save(common.genes, file ="topHVGNBNeuHB.RData")

q()
