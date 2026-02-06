## signal of SE based on sum of counts of constituent peaks
suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(ArchR)
})
addArchRThreads(threads = 8)

setwd("~/MBsnANA/HBana/ATACana/Outs/")


#
## load the subsetted peak matrix ------
## get the metadata and subset it
load("extrafiles/plotdataHB_ATACv6_0404255step2.RData")

cells <- readLines("extrafiles/CellforAUCscoresubsetbylv2.txt")
plot.data.ATAC <- plot.data.ATAC[cells,]

#load the project and peak matrix
proj <- loadArchRProject(path = "comATACx/")

##extract PeakMatrix for selected cells
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)

##robuset peaks to subset to
peaks <- read.delim("PeakCom/AllPeaks_robust.fil.sorted.bed", header = FALSE)
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
## also probably a good idea to normalize read counts
scale_factor <- colSums(PM)/10000

## --------
## obtain count matrix where each SE signal is sum of signal from constituent peaks

load("Superenhancer/rose/ALLSElist.RData")

allses <- metase$ID
mdt <- c()
for(x in allses){
  sel_cres <- all_se_cres[[x]]
  rs <- colSums(PM[sel_cres,])
  rs <- rs /scale_factor
  mdt <- rbind(mdt,rs)
  rm(sel_cres,rs)
}
row.names(mdt) <- allses
mdt[1:4,1:4]

save(mdt, plot.data.ATAC, file = "Superenhancer/SEcounts.RData")
## marker SEs per class ---
marks <- scran::findMarkers(mdt, groups=plot.data.ATAC$Class,
                            direction="up",
                            test.type="wilcox",pval.type="some")

cls <- names(marks)
sigs <- c()
for(x in cls){
  del <- marks[[x]]
  del <- row.names(del)[1:5]
  sigs <- c(sigs,del)
  rm(del)
  
  
}
sigs <- unique(sigs)

annlv2 <- unique(plot.data.ATAC$Annotationlv2_step2)

msd <- c()
for(x in annlv2){
  cells <- row.names(plot.data.ATAC)[plot.data.ATAC$Annotationlv2_step2==x]
  mdtx <- rowMeans(mdt[,cells])
  msd <- rbind(msd,mdtx)
  rm(cells,mdtx)
  
}
row.names(msd) <- annlv2
msd[1:4,1:4]

ann <- read.delim("~/MBsnANA/Diemana/scepro/Classcluster.txt")
row.names(ann) <- ann$Cluster

colr <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
colx <- colr$Color 
#annlv2col <- colorRampPalette(pal_rickandmorty()(12))(157)
#names(annlv2col) <- ann$Cluster
names(colx) <- row.names(colr)

# ann_col <- list(Class=colx,
#                 Cluster=annlv2col)
ann_col <- list(Class=colx)
ann$Cluster <- NULL

ann2 <- metase
ann2$SE <- ann2$ID <- ann2$Len <- NULL

#range 0 to 1----
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}

msd2 <- apply(msd,2,range_01)
row.names(msd2) <- row.names(msd)

mc <- colorRampPalette(c("white","grey75","black"))(101)
mb <- (seq(1:101)-1)/100

hp <- pheatmap::pheatmap(msd2[,sigs], 
                         clustering_method = "ward.D2",
                         #clustering_distance_cols = "correlation",
                         #clustering_distance_rows = "correlation",
                         #annotation_row = ann,
                         #annotation_col = ann2,
                         #annotation_colors = ann_col,
                         color=mc, breaks=mb,
                         border_color = NA,
                         fontsize = 5)
pdf("SEmarkerv2.pdf", width = 15, height = 12, pointsize = 10)
print(hp)
dev.off()



 ##