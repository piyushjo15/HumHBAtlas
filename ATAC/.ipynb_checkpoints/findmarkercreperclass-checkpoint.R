## getting peaks called on cluster_subset
suppressPackageStartupMessages({
  library(ArchR)
  library(scran)
  library(scater)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(dplyr)
  library(tidyr)
})
addArchRThreads(8)
setwd("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/")

## part 1. load all the data ------------
### 1.1 load plotdata and subset to max 200 cells
load("plotdataHBATAC_210725.RData")
proj <- loadArchRProject(path = "~/MBsnANA/HBana/ATACana/Outs/comATACx/")

cells1 <- getCellNames(proj)
cells2 <- row.names(plot.data)
cells1 <- cells1[cells1 %in% cells2]
plot.data <- plot.data[cells1,]
proj <- proj[cells1,]
table(row.names(plot.data)==getCellNames(proj))

proj$Ann <-plot.data$Class
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Ann",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 5000,k = 300
)
save(markersPeaks, file = "ATAC/markerpeaks_Class.RData")
q()
setwd("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/")
load("ATAC/markerpeaks_Class.RData")

markerList <- getMarkers(markersPeaks, cutOff = "FDR < 0.01 & (Log2FC >=2 | AUC >= 0.52)")
clusters <- names(markerList)
npeaks <- lengths(markerList)
all_peaks <- read.delim("~/MBsnANA/HBana/ATACana/Outs/PeakCom/AllPeaks_robust.fil.sorted.bed",header = FALSE,row.names = 4)
head(all_peaks)

marker_peaks <- c()
##focusing on max top 10k per class
for(x in clusters){
  del <- data.frame(markerList[[x]])
  del <- del[,c(1,3,4)]
  colnames(del) <- c("chr","start","end")
  id <- paste0(del$chr,":",del$start,"-",del$end)
  del$peak <- id
  keep <- del$peak %in% row.names(all_peaks)
  del <- del[keep,]
  del$Class <- x
  if(npeaks[x]>10000){
    del <- del[1:10000,]
  }
  marker_peaks <- rbind(marker_peaks,del)
  rm(del)
}
rm(x)
table(marker_peaks$Class)
save(marker_peaks, file="ATAC/markerpeaks_class_com.RData")

### proportion of peaks per chromosomes

## want to find the proportion of peaks that lie in a chromsome per class
df <- t(as.matrix(table(marker_peaks$Class,marker_peaks$chr)))
df[1:4,1:4]

rs <- rowSums(df)
cs <- colSums(df)
df_row_norm <- df / rs
df_col_norm <- t(t(df) / cs)
df_col_norm <- df_col_norm *100

df_row_norm[1:4,1:4]
df_col_norm[1:4,1:4]

library(pheatmap)

mc <- colorRampPalette(c("white","navajowhite","orange1","tan4"))(50)
mb <- seq(0, max(df_col_norm), length.out=50)
chr <- paste0("chr",seq(1:22))

hp <- pheatmap(df_col_norm[chr,],
         cluster_method="Ward.D2",
         cluster_rows=FALSE,
         color=mc,breaks=mb)

#pdf("Figures/ClassCREsperchromosome.pdf", width = 4, height = 5, pointsize = 20)
#print(hp)
#dev.off()

## significance test, using binomial
mat <- df

total_events <- sum(mat)
col_prop <- colSums(mat) / total_events
row_tot <- rowSums(mat)

n <- nrow(mat)
m <- ncol(mat)

pvals_binom <- matrix(NA, n, m)

for (i in 1:n) {
  for (j in 1:m) {
    x  <- mat[i, j]        # observed
    N  <- row_tot[i]       # total events in that class
    p  <- col_prop[j]      # expected proportion
    
    # one-sided enrichment test
    pvals_binom[i, j] <- binom.test(x, N, p = p, alternative = "greater")$p.value
  }
}

#fdr
pvals_fdr <- matrix(
  p.adjust(as.vector(pvals_binom), method = "fdr"),
  nrow = nrow(pvals_binom),
  ncol = ncol(pvals_binom)
)

pval_sig <- matrix(0,nrow(pvals_fdr),ncol(pvals_fdr))
pval_sig[pvals_fdr<0.01] <- 1

row.names(pval_sig) <- row.names(df)
colnames(pval_sig) <- colnames(df)


#pdf("Figures/ClassCREsperchromosomeFDRsig.pdf", width = 3.5, height = 5, pointsize = 20)
#print(hp2)
#dev.off()

##ggplot bubbple plot
df2 <-reshape2::melt(df_col_norm[chr,])
df3 <-reshape2::melt(pval_sig[chr,])
df2$Var1 <- as.character(df2$Var1)
df2$Var2 <- as.character(df2$Var2)
df3$Var1 <- as.character(df3$Var1)
df3$Var2 <- as.character(df3$Var2)

row.names(df2) <- paste0(df2$Var1,"_",df2$Var2)
row.names(df3) <- paste0(df3$Var1,"_",df3$Var2)
colnames(df2) <- c("chr","Class","Prop")
table(row.names(df2)==row.names(df3))
df2$Sig <- ifelse(df3$value==1,TRUE,FALSE)


df2$chr <- factor(df2$chr,levels=rev(chr))
cls<- c("Progenitors","Neuroblasts_I","Neuroblasts_II","Neuroblasts_III",
            "GCUBC","GC","CBINT","Purkinje", "Neurons_I", "Neurons_II","Glioblasts",
            "PreO","MOG","MixedGlia","MES_VAS","Immune")
df2$Class <- factor(df2$Class,levels=cls)

p <- ggplot(df2, aes(x = Class, y = chr)) +
  geom_point(aes(size = Prop, fill = Prop), shape = 21, color = "white") +
  geom_point(
    data = subset(df2, Sig),
    aes(size = Prop),
    shape = 21, fill = NA, color = "black", stroke = 1.1
  ) +
  scale_size_continuous(range = c(0,12), limits = c(0,max(df2$Prop))) +
  scale_fill_gradientn(colors = mc, limits = c(0,max(df2$Prop))) +
  theme_minimal()+
theme( axis.text.x = element_text(colour = "black",size=10,angle = 45,hjust=1),
     axis.text.y = element_text(colour = "black",size=10))
pdf("Figures/ClassCREsperchromosome.pdf", width = 6, height = 6, pointsize = 10)
print(p)
dev.off()
q()

DIR="~/MBsnANA/HBana/ATACana/Outs/markerCREs/Class/"
cls <- unique(marker_peaks$Class)
for(x in cls){
    del <- marker_peaks[marker_peaks$Class==x,]
    del <- del[1:2000,]
    write.table(del, file=paste0(DIR,x,"_CREs.bed"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    rm(del)
    }

## obtaining hg19 version of the marker peaks for GWAS

all_peaks_hg19 <- read.delim("~/MBsnANA/HBana/ATACana/Outs/PeakCom/AllPeaks_robust.fil.hg19v2.bed",header = FALSE)
rn <- row.names(all_peaks_hg19) <- all_peaks_hg19$V4
head(all_peaks_hg19)


for(x in clusters){
  del <- data.frame(markerList[[x]])
  del <- del[,c(1,3,4)]
  colnames(del) <- c("chr","start","end")
  id <- paste0(del$chr,":",del$start,"-",del$end)
  del$peak <- id
  keep <- del$peak %in% row.names(all_peaks)
  del <- del[keep,]
  del$Class <- x
  if(npeaks[x]>20000){
    del <- del[1:20000,]
  }
  id <- del$peak
  ##hg19
  rnx <- rn[rn %in% id]
  write.table(all_peaks_hg19[rnx,],paste0("~/MBsnANA/HBana/ATACana/Outs/GWAS/markersv4/",x,"_CREs.bed"),
              sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

  rm(del)
}
rm(x)

##
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR < 0.01 & (Log2FC >=2 | AUC >= 0.52)",
  transpose = TRUE
)
mdt <- heatmapPeaks@matrix
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#add background peaks----
proj <- loadArchRProject(path = "comATACx/")
proj <- addBgdPeaks(proj, force = T)

load("comATACx/PeakSet.RData")
head(PeakSet)
table(PeakSet$Reproducibility)
