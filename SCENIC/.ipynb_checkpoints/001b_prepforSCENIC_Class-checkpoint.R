## Preparing a loom file for each Class
suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(Matrix)
  library(dplyr)
  library(tidyr)

})

#Directories
DIR_RNA <-"~/MBsnANA/Diemana/scepro/"
DIR_ATAC <-"~/MBsnANA/HBana/ATACana/Outs/"
DIR_OUT <-"~/MBsnANA/HBana/SCENICana/SCENICOut/Class"
setwd(DIR_OUT)

args <- commandArgs(trailingOnly = TRUE)
print(paste0("Generating data for loom file pyscenic run for ..", args))

## 1. --------
#loading counts file
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
#load hvg data
load(paste0(args,"/",args,"_HVGnTF.RData"))
if(length(hvg)>2000){
  hvg <- hvg[1:2000]
  
}
hvg <- unique(c(hvg,heg,tfs_sel))

## obtain plotdata
load(paste0(DIR_RNA,"NMF_MP/plotdataHBforNMFMP.RData"))
plot.data <- plot.data[plot.data$Annotationlv1==args,]
com.sce <- com.sce[,row.names(plot.data)]

#doing for each splits -----
splts <- c("A","B")
for(x in splts){
  print(paste0("Running for split:",x))
  pds <- plot.data[plot.data$SplitNMF==x,]
  # downsampling
  tb <- data.frame(table(pds$Annotationlv2))
  cl <- as.character(tb$Var1)[tb$Freq>3000]
  keep <- pds$Annotationlv2 %in% cl
  table(keep)
  pds1 <- pds[keep, ]
  pds2 <- pds[!keep, ]
  pds3 <- c()
  rm(keep,tb)
  ##downsampling large annotations to max 3k cells
  for(y in cl){
    del <- pds1[pds1$Annotationlv2==y,]
    rn <- row.names(del)
    del2 <- del[sample(rn,3000),]
    pds3 <- rbind(pds3,del2)
    rm(del,rn,del2)
  }
  rm(y)
  pds4 <- rbind(pds2,pds3)
  
  ##loom file
  print("Creating data for loom file")
  exprMat <- round(counts(com.sce)[hvg,row.names(pds4)])
  cellInfo <- pds4
  
 
  ##save this temporarly to build loom
  save(exprMat,cellInfo, file = paste0(args,"/",args,"_",x,"_4loom.RData"))
  
  rm(exprMat,cellInfo,pds,pds1,pds2,pds3,pds4,cl)
 
  
}

##created loom using 001c_cnvrt2loom.R

q()
## old for differen t number of splits
