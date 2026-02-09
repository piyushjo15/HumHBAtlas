##finding differentially accessible genes in OPC compartment of individual tumor sample
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(tidyverse)
})

##Directories 
DIR_IN="~/MBsnANA/HGGana/Mija_DIPG/"
DIR_NOR <- "~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/"
sam <- commandArgs(trailingOnly = TRUE)
setwd(DIR_IN)

## lod HVG calculated for tumor
load("DMG_OPC_HVG.RData")
load("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/Nor_OPC_HVG.RData")
genes <- unique(c(hvg_tum,hvg_nor))


## GLM GAMP poi-----------------
library(glmGamPoi)
## this is compiled log transformed expression data for the entire hindbrain
load("combinedcountsMijaDIPG10x.RData")
mdt <- round(counts(com.sce))
rm(com.sce)
##load metadata, combined
load("plotdataMijaDIPG10x_ANN.RData")
table(row.names(plot.data)==colnames(mdt))
## extract Oligo cells from sample
plot.data <- plot.data[plot.data$Batch==sam,]

## extract Oligo cells

cl <- grep("Oligo",plot.data$ANN2)
pd_tum <- plot.data[cl,]
mdt_tum <- mdt[genes,row.names(pd_tum)]
pd_tum$Class <- "Tumor"
pd_tum$ANN <- pd_tum$ANN2
pd_tum$ANN2 <- pd_tum$ANN3

rm(mdt,plot.data)
## normal
load(paste0(DIR_NOR,"plotdataHBRNA_241125.RData"))

# subset to glial lineage
cl <-  c("PreO:OPC_PrN" ,"PreO:OPC_PNAD" ,"PreO:NFOL","PreO:COP")
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

load(paste0(DIR_NOR,"combinedcountsHB.RData"))
com.sce <- com.sce[,row.names(pd_nor)]
mdt_nor <- round(counts(com.sce[genes,]))
pd_nor$Class <- "Normal"
pd_nor$Batch <- "Normal"
pd_nor$ANN <- pd_nor$Annotationlv2
pd_nor$ANN2 <- pd_nor$Annotationlv3

rm(com.sce,plot.data)


## merge
mdt <- cbind(mdt_nor,mdt_tum)
cls <- c("Batch","Class","ANN","ANN2")
pd <- rbind(pd_nor[,cls],pd_tum[,cls])
ann <- pd$Class
setwd(paste0("/home/p541i/MBsnANA/Diemana/scepro/RECORDR/",sam))
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
save(res, file="degOPCvsTumOPC_glm_new.RData")
q()
