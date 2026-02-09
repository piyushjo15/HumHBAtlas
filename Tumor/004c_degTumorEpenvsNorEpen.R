suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(tidyverse)
    library(glmGamPoi)

})

##Directories 
DIR_IN="~/PFA/"
sam <- commandArgs(trailingOnly = TRUE)
setwd(DIR_IN)

## lod HVG calculated for tumor
load("PFA_Epen_HVG.RData")
load("~/RNA/Nor_Epen_HVG.RData")
genes <- unique(c(hvg_tum,hvg_nor))


##1. load required data -------

## normal
load("~/RNA/plotdataHBRNAv2.RData")

# subset to glial lineage
cl <- c("Glioblasts:Ependoblast","MixedGlia:Ependymal")

pd_nor <- plot.data[plot.data$Annotationlv2 %in% cl,]
## sample down to  3 k cells
cl <- unique(pd_nor$Annotationlv2)
cells <- c()
for(x in cl ){
  keep <- row.names(pd_nor)[pd_nor$Annotationlv2==x]
  if(length(keep)>3000){
    keep <- sample(keep,3000)
  }
  cells <- c(cells,keep)
  rm(keep)
}
rm(x)
pd_nor <- pd_nor[row.names(pd_nor) %in% cells,]
cells <- row.names(pd_nor)

load("~/RNA/combinedcountsHB.RData"))
com.sce <- com.sce[,row.names(pd_nor)]
mdt_nor <- round(counts(com.sce))
genes_n <- row.names(mdt_nor)

pd_nor$Class <- "Normal"
pd_nor$Batch <- "Normal"
pd_nor$ANN <- pd_nor$Annotationlv2
pd_nor$ANN2 <- pd_nor$Annotationlv3

rm(com.sce,plot.data)


## load tumor expression data
## this is compiled log transformed expression data for the entire hindbrain
load("combinedcountsGJPFA.RData")
mdt <- round(counts(com.sce))
rm(com.sce)

genes_t <- row.names(mdt)

genesx <- intersect(genes_t,genes_n)
genes <- intersect(genes,genesx)

mdt_nor <-mdt_nor[genes,]

##load metadata, combined
load("plotdataGJPFA_ANN.RData")
table(row.names(plot.data)==colnames(mdt))

## extract Ependymal cells
cl <- grep("Ependymal",plot.data$ANN2)
pd_tum <- plot.data[cl,]
mdt_tum <- mdt[genes,row.names(pd_tum)]
pd_tum$Class <- "Tumor"
pd_tum$ANN <- pd_tum$ANN2
pd_tum$ANN2 <- pd_tum$ANN3

rm(mdt,plot.data)


## merge
mdt <- cbind(mdt_nor,mdt_tum)
cls <- c("Batch","Class","ANN","ANN2")
pd <- rbind(pd_nor[,cls],pd_tum[,cls])
ann <- pd$Class
setwd(paste0("/home/PFA/",sam))
mdt <- as.matrix(mdt)
fit <- glm_gp(
  mdt,
  design = ~ ann,
  col_data = pd
)
res <- test_de(
  fit,
  contrast = "annTumor"
)
row.names(res) <- res$name
save(res, file="degEpenvsTumEpen_glm_new.RData")
q()
