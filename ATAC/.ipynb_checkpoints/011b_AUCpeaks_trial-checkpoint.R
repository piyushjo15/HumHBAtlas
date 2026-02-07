##This scripts calculates AUC score of GRN set for specific lineage 
suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(AUCell)
  library(BiocParallel)
  library(ArchR)
})

addArchRThreads(threads = 8)

args <- commandArgs(trailingOnly = TRUE)
print(paste0("Calcualting AUC for class: ",args))
##
DIR <- "~/ATACana/Outs/"
setwd(DIR)

## 1. load peak associated with metaprograms ---------
DIR_INP <- "~//FetalBrain/markerpeaks/"
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
load("extrafiles/plotdataHB_ATACv6_0404255step2.RData")
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$Class==args,]
cells <- row.names(plot.data.ATAC)
## 3. load peak Matrix-----
proj <- loadArchRProject(path = "comATACx_4peaks/")

##extract PeakMatrix for selected cells
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)
##robuset peaks to subset to
peaks <- read.delim("~/MBsnANA/HBana/ISM/Outs/FetalBrain/All_peak_FetalBrain.500bp.strtmin1.bed", header = FALSE)
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
save(plot.data.ATAC,scr, file =paste0("CREs/pdAUC_FB_CRE_",args,".RData"))

q()
### analysis for HB Topic CRE ------
cls <- readLines("~/MBsnANA/HBana/ATACana/HBATAC_Scripts/RNA_Class.txt")

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

load("CREs/MergedpdAUC_FB_CRE.RData")
df <- data.frame(table(pd$Annotationlv2_step2))
df <- as.character(df$Var1)[df$Freq>25]
keep <- pd$Annotationlv2_step2 %in% df
pd <- pd[keep,]
mdt <- mdt[keep,]
annlv2 <- unique(pd$Annotationlv2_step2)

msd <- c()
for(x in annlv2){
  cells <- row.names(pd)[pd$Annotationlv2_step2==x]
  mdtx <- colMeans(mdt[cells,])
  msd <- rbind(msd,mdtx)
  rm(cells,mdtx)
  
}
row.names(msd) <- annlv2
msd[1:4,1:4]

msd2 <- apply(msd,2,scale)
row.names(msd2) <- row.names(msd)

msd2[msd2>2] <- 2
msd2[msd2<(-2)] <- (-2)


mc <- colorRampPalette(c("deepskyblue4","beige","magenta3"))(101)
mb <- c(seq(-2, 0, length.out=ceiling(100/2) + 1),
              seq(2/100, 2, length.out=floor(100/2)))


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

sel_mp <- c("MP35","MP9","MP6","MP51","MP53","MP54")
sel_mp <- readLines()

ct_ord <-readLines("~/MBsnANA/HBana/ISM/Outs/FetalBrain/Celltype_order.txt")
mp_ord <-readLines("~/MBsnANA/HBana/ISM/Outs/FetalBrain/MP_order.txt")

mdt_mp <- mdt_mp[,mp_ord]
corx <- cor(mdt_fb,mdt_mp)

mc <- colorRampPalette(c("blue4","blue2","white","orangered2","orangered3"))(201)
mb <- ((seq(1:201)-1)/100)-1


cl <- c("Neur_Purk_4","Neur_Purk_2" ,"Neur_Purk_3","Neur_Purk_5","Neur_Purk_1","Neur_Purk_6",
        "Mic_active_1","PVM", "Mic_2","Mic_3","Mic_dividing","Mic_5","Mic_6","Mic_active_2","Mic_1",
        "Immune_uncertain","Mic_4","Schwl_2", "Schwl_3","PreOPC_1","Schwl_1",
        "T_Cell_immature","Neur_CB_GNP_IPC_2","Neur_CB_GNP_IPC_4","Neur_CB_GNP_IPC_1",
        "Neur_CB_GNP_IPC_3","Neur_Purk_NBL_1","Neur_Purk_NBL_2","Neur_Purk_prog",
        "COP","OPC_6","OPC_2","OPC_1","OPC_3","PreOPC_3","PreOPC_4","OPC_4","OPC_5",
        "PreOPC_2")
hp <- pheatmap::pheatmap(corx[ct_ord,], 
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         clustering_method = "ward.D2",
                         color=mc, breaks=mb,
                         legend = FALSE,
                         border_color = "grey",
                         fontsize = 12)
png("Figures/FetalbrainCorMPlabs.png",units = "in",width = 6,height = 8,res = 300)
print(hp)
dev.off()
####-----

ann <- read.delim("~/MBsnANA/Diemana/scepro/annlv2_x.txt")
row.names(ann) <- ann$Annlv2
colr <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
colx <- colr$Color 
annlv2col <- colorRampPalette(pal_rickandmorty()(12))(157)
names(annlv2col) <- ann$Annlv2
names(colx) <- row.names(colr)
ann_col <- list(NewClass=colx,
                Annlv2=annlv2col)

hp <- pheatmap::pheatmap(msd2, 
                         clustering_method = "average",
                         clustering_distance_cols = "correlation",
                         clustering_distance_rows = "correlation",
                         annotation_row = ann,
                         annotation_colors = ann_col,
                         color=mc, breaks=mb,
                         border_color = NA,
                         fontsize = 8)

pdf("MPbyannlv2enr.pdf", width = 16, height = 16, pointsize = 10)
print(hp)
dev.off()

##quantile normalize
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
msd2 <- apply(msd,2,range_01)

msd3 <- preprocessCore::normalize.quantiles(msd2)
colnames(msd3) <- colnames(msd)
row.names(msd3) <- row.names(msd)

mc <- colorRampPalette(c("white","grey75","black"))(101)
mb <- (seq(1:101)-1)/100


hp <- pheatmap::pheatmap(msd3, 
                         clustering_method = "ward.D2",
                         annotation_row = ann,
                         annotation_colors = ann_col,
                         color=mc, breaks=mb,
                         border_color = NA,
                         fontsize = 8)

pdf("FBmarkerannlv2enrQN.pdf", width = 16, height = 16, pointsize = 10)
print(hp)
dev.off()
save(msd3, file = "CREs/FB_msd3.RData")
##