## This script compares AUC scores obtained for cell cluster specific 
## marker peaks obtained in Fetal brain atlas (Mannens Nature 2024)
## to best matching meta-regulatory programs

## first to calculate AUC values per marker-peak set, the code was run per class
## else it needed too much memory. 
## Then the AUC values were combined to obtain integrated values
## Then a correlation was obtaine dbetween AUC values for markepeaks and rMPs
suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(AUCell)
  library(BiocParallel)
  library(ArchR)
})

addArchRThreads(threads = 8)

Cls <- commandArgs(trailingOnly = TRUE)
print(paste0("Calcualting AUC for class: ",Cls))
##
DIR <- "~/ATACana/Outs/"
setwd(DIR)

## 1. load peak associated with Fetal brain atlas ---------
DIR_INP <- "~/DeepHB/Outs/FetalBrain/markerpeaks/"
#for markers from Fetal brain
mps <- list.files(DIR_INP)
mps <- gsub("_peaks.bed","",mps, fixed = TRUE)
mp_list <- list()
for(x in mps){
 del <- read.delim(paste0(DIR_INP,x,"_peaks.bed"),header = FALSE)
 mp_list[[x]] <- del$V4
 rm(del)
}



## 2. load metadata -------
##here I am not subsetting or selecting for any annotation, just for class
load("plotdataHB_ATACv6.RData")
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$Class==Cls,]
cells <- row.names(plot.data.ATAC)
## 3. load peak Matrix-----
## a separate ArchRproject was created to add peaks identified in Fetal brain atlas
## to obtain Peak Matrix for it
proj <- loadArchRProject(path = "comATACx_4Fetpeaks/") 

##extract PeakMatrix for selected cells
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)
##robuset peaks to subset to
## The peaks were adujust by one base 
peaks <- read.delim("~/DeepHB/Outs/FetalBrain/All_peak_FetalBrain.500bp.strtmin1.bed", header = FALSE)
peaks$V2 <- peaks$V2+1
peaks <- paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)

region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)

##get Matrix
PM <- assay(PeakMatrix)
row.names(PM) <- row.names(region_names)
PM[1:4,1:4]
PM <- PM[peaks,cells] ###ver important to make sure cell order is same
dim(PM)
rm(PeakMatrix,peaks)


###AUCell-----
ms <- AUCell_run(PM, mp_list, aucMaxRank=nrow(PM)*0.10, normAUC = FALSE,
                 BPPARAM=BiocParallel::MulticoreParam(9))
scr <- t(assay(ms))
save(plot.data.ATAC,scr, file =paste0("CREs/pdAUC_FB_CRE_",Cls,".RData"))

q()
### analysis for HB Topic CRE ------
cls <- readLines("~/RNA/RNA_Class.txt")

pd <- c()
mdt <- c()
for(x in cls){
  load(paste0("CREs/pdAUC_FB_CRE_",x,".RData"))
  mdt <- rbind(mdt,scr)
  pd <- rbind(pd,plot.data.ATAC)
  rm(plot.data.ATAC,scr)
  
}
table(row.names(pd) %in% row.names(mdt))
mdt <- mdt[row.names(pd),]
save(mdt, pd, file = "CREs/MergedpdAUC_FB_CRE.RData")


## correlation of Fetal brain with MP cre---
load("CREs/MergedpdAUC_FB_CRE.RData")
mdt_fb <- mdt
load("CREs/MergedpdAUC_MPCRE.RData")
mdt_mp <- mdt2
dim(mdt_fb)
dim(mdt_mp)
rm(mdt,mdt2)
rn <- intersect(row.names(mdt_fb),row.names(mdt_mp))
mdt_fb <- mdt_fb[rn,]
mdt_mp <- mdt_mp[rn,]

mp_ord <- c("rMP33","rMP53","rMP54","rMP52","rMP55",
            "rMP24","rMP27","rMP32","rMP20","rMP37","rMP39")

ct_ord <- c("Neur_Purk_4","Neur_Purk_2" ,"Neur_Purk_3","Neur_Purk_5","Neur_Purk_1","Neur_Purk_6",
        "Mic_active_1","PVM", "Mic_2","Mic_3","Mic_dividing","Mic_5","Mic_6","Mic_active_2","Mic_1",
        "Immune_uncertain","Mic_4","Schwl_2", "Schwl_3","PreOPC_1","Schwl_1",
        "T_Cell_immature","Neur_CB_GNP_IPC_2","Neur_CB_GNP_IPC_4","Neur_CB_GNP_IPC_1",
        "Neur_CB_GNP_IPC_3","Neur_Purk_NBL_1","Neur_Purk_NBL_2","Neur_Purk_prog",
        "COP","OPC_6","OPC_2","OPC_1","OPC_3","PreOPC_3","PreOPC_4","OPC_4","OPC_5",
        "PreOPC_2")
mdt_mp <- mdt_mp[,mp_ord]
corx <- cor(mdt_fb,mdt_mp)

mc <- colorRampPalette(c("blue4","blue2","white","orangered2","orangered3"))(201)
mb <- ((seq(1:201)-1)/100)-1



hp <- pheatmap::pheatmap(corx[ct_ord,], 
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         clustering_method = "ward.D2",
                         color=mc, breaks=mb,
                         legend = FALSE,
                         border_color = "grey",
                         fontsize = 12)
