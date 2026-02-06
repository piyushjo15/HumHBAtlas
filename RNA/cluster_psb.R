suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
})

#library(scuttle)
#args <- "Neurons"
setwd("~/MBsnANA/Diemana/scepro/")

# ## psb ------------
load("combinedcountsHB.RData")
# 
load("plotdataHBTempANN250305.RData")
plot.data <- plot.data[plot.data$Annotationlv3!="ND",]
#tb <- data.frame(table(plot.data$Annotationlv3))
#tb <- as.character(tb$Var1)[tb$Freq>24]
#keep <- plot.data$Annotationlv3 %in% tb
#plot.data <- plot.data[keep,]

## pseudobulk by class----------
## subset
cells <- readLines("CellformodscoreHBv3.txt")
plot.data <- plot.data[row.names(plot.data) %in% cells,]
del <- com.sce[,row.names(plot.data)]
del <- aggregateAcrossCells(del,id=plot.data$Annotationlv1)
del <- logNormCounts(del)
save(del, file ="HBAnnlv1psb_280425.RData")
rm(del)
# ## get HVG by Class block
# load("lgcountsHB.RData")
# mdt <- mdt[,row.names(plot.data)]
# dec <- modelGeneVar(mdt,block=plot.data$Annotationlv1)
# topHVG <- getTopHVGs(dec)
# load("~/MBsnANA/Diemana/scepro/genesnoMTRibo.RData")
# topHVG <- topHVG[topHVG %in% genes]
# load("~/MBsnANA/Diemana/scepro/HGNCsexchr.RData")
# topHVG <- topHVG[!(topHVG %in% SX)]
# save(dec, topHVG, file = "HBClassHVGblckbyclass.RData")
# #q()
# 
# 
# # df3 <- data.frame(table(plot.data$Annotationlv2,plot.data$Class))
# # df3 <- df3[df3$Freq>0,]
# # write.table(df3, file = "Annlv2Class.txt",quote = FALSE,row.names = FALSE,
# #             sep = "\t")
## lv2---------
# cells <- c()
# cls <- unique(plot.data$Annotationlv2)
# for(x in cls){
#   cells_sel <- row.names(plot.data)[plot.data$Annotationlv2==x]
#   if(length(cells_sel)>2000){
#     cells_sel <- sample(cells_sel,2000)
#     cells <- c(cells,cells_sel)
#   } else {
#     cells <- c(cells,cells_sel)
#   }
#   rm(cells_sel)
# }
# rm(x)
# pd <- plot.data[cells,]
del <- com.sce[,row.names(plot.data)]
del <- aggregateAcrossCells(del,id=plot.data$Annotationlv2)
del <- logNormCounts(del)
save(del, file ="HBAnnlv2psb_280425.RData")
q()
##subsample and sum approach at lv3--------
cells <- c()
cls <- unique(plot.data$Annotationlv3)
for(x in cls){
  cells_sel <- row.names(plot.data)[plot.data$Annotationlv3==x]
  if(length(cells_sel)>1000){
    cells_sel <- sample(cells_sel,1000)
    cells <- c(cells,cells_sel)
  } else {
    cells <- c(cells,cells_sel)
  }
  rm(cells_sel)
}
rm(x)
pd <- plot.data[cells,]
del <- com.sce[,row.names(pd)]
del <- aggregateAcrossCells(del,id=pd$Annotationlv3)
del <- logNormCounts(del)
save(del, file ="HBAnnlv3psb_sub_del.RData")
q()

#mean
del <- counts(com.sce[,row.names(plot.data)])
rm(com.sce)
cls <- unique(plot.data$Annotationlv3)
del2 <- c()
for(x in cls){
  keep <- plot.data$Annotationlv3==x
  del3 <- apply(del[,keep],1,median)
  del2 <- cbind(del2,del3)
  rm(keep,del3)
}
row.names(del2) <- row.names(del)
colnames(del2) <- cls
del <- SingleCellExperiment(list(counts=del2))
del <- logNormCounts(del)
rm(del2)
save(del, file ="HBAnnlv3psb_median.RData")
rm(del)
# #sum
# del <- com.sce[,row.names(plot.data)]
# del <- aggregateAcrossCells(del,id=plot.data$Annotationlv3)
# del <- logNormCounts(del)
# save(del, file ="HBAnnlv3psb.RData")
# q()

load("lgcountsHB.RData")
mdt <- mdt[,row.names(plot.data)]
del <- c()
for(x in cls){
  keep <- plot.data$Annotationlv3==x
  del2 <- apply(mdt[,keep],1,median)
  del <- cbind(del,del2)
  rm(keep,del2)
}
row.names(del) <- row.names(mdt)
colnames(del) <- cls
save(del,file ="HBAnnlv3psb_median_lg.RData")
q()
library(ggsci)
#
## Trying consensus clustering to find number of classes ----------
load("~/MBATAC/MBnewana/HGNCsexchr.RData")
load("HBAnnlv3psb.RData")
#com.hvg <- hvg[!(hvg %in% SX)]
load("topHVGcomhvgnoRiboMTSX.RData")
psb <- logcounts(del)
psb <- psb[topHVG[1:4000],]
dim(psb)
cl <- c("Neurons_3:Dopamine_Cer:x",
        "Neurons_3:LHX2_PRDM8:x")
clx <- colnames(psb) %in% cl
psb2 <- t(apply(psb[,!clx], 1, scale))
colnames(psb2) <- colnames(psb)[!clx]
results = ConsensusClusterPlus(psb2,maxK=20,reps=100,pItem=0.8,pFeature=0.75,
                               innerLinkage="complete",finalLinkage="complete",
                               title="CCHCWD24kS",clusterAlg="hc",
                               distance="pearson",seed=123,plot="png")

cluster = data.frame(X=results[[18]][["consensusClass"]])
cluster$X <- paste0("CL",cluster$X)
write.table(cluster, file = "CCHCom4kS_CO18.txt",sep = "\t", quote = FALSE)


icl = calcICL(results,title="CCHCWD24kS",plot="png")

icl[["itemConsensus"]][1:5,]
save(results, icl, file = "CCHCWD24kS.RData")

### plotting based on CC tree--------
load("HBAnnlv3psb.RData")
load("CCHCWD24kS.RData")
psb <- logcounts(del)
psb <- psb[cand,]

ann <- read.delim("Classcluster.txt")
row.names(ann) <- ann$Ann3
ann$Ann3 <- ann$Ann2 <- NULL

col1 <- colorRampPalette(pal_cosmic()(10))(length(unique(ann$Ann1)))
col2 <- colorRampPalette(pal_simpsons()(14))(length(unique(ann$Class)))
names(col1) <- unique(ann$Ann1)
names(col2) <- unique(ann$Class)
ann_col <- list(Ann1=col1,
                Class=col2)

psb2 <- t(apply(psb, 1, zscore_c))

psb2[psb2>2.5]<- 2.5
psb2[psb2<(-2.5)]<- (-2.5)

mc <- colorRampPalette(c("deepskyblue4","beige","magenta3"))(101)
mb <- c(seq(-2.5, 0, length.out=ceiling(100/2) + 1),
        seq(2.5/100, 2.5, length.out=floor(100/2)))

## this plots the tree from CC
X=results[[18]]
hc=X$consensusTree
del2 <- colnames(psb2)[hc$order]
write(del2, file = "orderofclustersinCC.txt")
library(grid)
hp <-pheatmap(psb2, 
              cluster_cols =  as.hclust(hc),
              color=mc, breaks=mb,
              border_color = NA,fontsize = 8, 
              show_rownames = FALSE,show_colnames = FALSE,
              annotation_col = ann, annotation_colors = ann_col)
tiff("CCHBpsbANNlv3.tiff", units="in", width=14, height=8, res=150)
print(hp)
dev.off()


##reorder based on given order

##rotated tree does not look good here
sel <- readLines("Annlv3ord.txt")
kk <- colnames(psb)
sel <- sel[(sel %in% kk)]
tcol <- hp$tree_col
tcoln <- dendextend::rotate(tcol,sel)

par(mfrow=c(2,1))
plot(trow)
tiff("HbTreeANNlv3.tiff", units="in", width=14, height=4, res=300)
plot(tcoln, cex=0.1)
dev.off()
##
ann <- read.delim("Classclusterv2.txt")
row.names(ann) <- ann$Annlv3
ann <- ann[,c("Annlv3","AnnLv2new","Class")]
ann <- ann[sel,]

anny <- read.delim("Classcol.txt",row.names = 1)
clscl <-anny$Color
names(clscl) <- row.names(anny)

d <- length(unique(ann$AnnLv2new))
lv2col <- colorRampPalette(pal_simpsons()(16))(d)
names(lv2col) <- unique(ann$AnnLv2new)

d <- length(unique(ann$Annlv3))
lv3col <- colorRampPalette(pal_igv()(50))(d)
names(lv3col) <- unique(ann$Annlv3)


ann_cols <- list(Class=clscl,
                 AnnLv2new=lv2col,
                 Annlv3=lv3col)


tiff("HbTreeANNlv3CM.tiff", units="in", width=14, height=6, res=300)
hp <- pheatmap::pheatmap(psb2[cand,sel], clustering_method = "ward.D2",
                         cluster_cols =  FALSE,
                         show_colnames = FALSE,
                         color=mc, breaks=mb,
                         border_color = NA,fontsize = 8, 
                         annotation_colors = ann_cols,
                         annotation_col = ann,
                         annotation_legend = FALSE)
dev.off()


##correlation--------
cl <- c("Neurons_3:Dopamine_Cer:x",
       "Neurons_3:LHX2_PRDM8:x")
clx <- colnames(psb) %in% cl
dmat <- dist(t(psb[,!clx]),method = "euclidean")
hc <- hclust(dmat, method = "ward.D2")
plot(hc,cex=0.5)

psb2 <- t(apply(psb, 1, scale))
colnames(psb2) <- colnames(psb)
cx <- cor(psb2,psb2, method = "pearson")
mc <- colorRampPalette(c("blue4","white","orangered4"))(101)
dmat <- as.dist(1 - cx) # Using 1 - correlation as distance
hc <- hclust(dmat, method = "complete")

cluster = data.frame(X=cutree(hc, k = 18))
cluster$X <- paste0("CL",cluster$X)
cluster$X <- paste0("CL",py$clusters)

write.table(cluster, file = "ann_tree_com_del.txt",sep = "\t", quote = FALSE)
#annotion
ann <- read.delim("Classcluster.txt")
row.names(ann) <- ann$ANN3
ann$ANN3 <- ann$ANN2 <- NULL

col1 <- colorRampPalette(pal_cosmic()(10))(length(unique(ann$ANN1)))
col2 <- colorRampPalette(pal_simpsons()(14))(length(unique(ann$Class)))
names(col1) <- unique(ann$ANN1)
names(col2) <- unique(ann$Class)
ann_col <- list(ANN1=col1,
                Class=col2)
hp <- pheatmap::pheatmap(psb2, 
                         clustering_method = "ward.D2",
                         clustering_distance_cols =  "correlation",
                         color=mc, breaks=mb,
                         border_color = NA,fontsize = 8, 
                         #cluster_rows=FALSE, 
                         show_rownames = FALSE, 
                         show_colnames = FALSE,
                         annotation_col = ann, annotation_colors = ann_col)
tiff("HBpsbANNlv3.tiff", units="in", width=14, height=12, res=150)
print(hp)
dev.off()



## this plots tree from HC but re-order columns as per order give,

tcol <- hp$tree_col
hp <-pheatmap(psb[cand,cl], scale = "row",
              cluster_cols =  as.hclust(tcol),
              color=mc, breaks=mb,
              border_color = NA,fontsize = 8, 
              #cluster_rows=FALSE, 
              #show_rownames = FALSE, 
              show_colnames = FALSE,
              legend = FALSE,annotation_legend = FALSE,
              annotation_col = ann, annotation_colors = ann_col)
tiff("HBpsbANNlv3_neuronalmarker.tiff", units="in", width=14, height=6, res=150)
print(hp)
dev.off()
# ##clustering by spearman clustering --------
# ##spearman clusters don't look good
# corx <- cor(psb2,psb2,method="spearman")
# dm=as.dist(1 - corx)
# hc <- hclust(dm,method = "ward.D2")
# 
# cluster = data.frame(X=cutree(hc, k = 18))
# cluster$X <- paste0("CL",cluster$X)
# write.table(cluster, file = "delhc.txt",sep = "\t", quote = FALSE)
# del3 <- hc$labels
# del3 <- del3[hc$order]
# write(del3, file = "delorderofclustersin.txt")

#markers
cand <- c("VIM","SOX2","NES",
          "MKI67","TOP2A","HES5","HES1",
          "OLIG1","OLIG2","PDGFRA",
          "ASCL1", "NEU4",
         "GPR17","SOX4","SOX11",
          "PDGFA","DOCK6","ITPR2",
          "MAP7","PLP1","CNP","MOG","MAG",
          "ELOVL7","RBFOX1","SVEP1","OPALIN",
          "GFAP","TNC","S100B","CD44","AQP1","EGFR",
          "LRRC4C","VAV3",
          "ADGRV1","CASC15","EPAS1","FABP7",
          "SLC1A3","SETBP1",
          "PTN","BCAN","CDH1","GLIS3",
          "FGFR3","AQP4","ELMO1",
          "CST3","SLC6A11","TRPM3","GJA1",
          "TNR", "PTGDS", "CRYAB", "SLC7A10",
          "AUTS2", "FREM2", "GRIA1",
          "CD24","SOX9","CFAP126","DNAI1",#EPENDYMAL
          "SLC17A7","SLC17A6","SLC17A8",
          "SLC6A2",
          "GAD1","GAD2",
          "SLC6A5","SLC32A1", "GABRB2","GABRG3",
          "TH","SLC18A2","DBH","SLC6A3","NR4A2","PPP1R1B",
          "TPH2","SLC6A4",
          "CHAT","SLC5A7","SLC18A3")
cand <- c(
          "SLC17A7","SLC17A6","SLC17A8",
          "SLC6A2",
          "GAD1","GAD2",
          "SLC6A5","SLC32A1", "GABRB2","GABRG3",
          "TH","SLC18A2","DBH","SLC6A3","NR4A2","PPP1R1B",
          "TPH2","SLC6A4",
          "CHAT","SLC5A7","SLC18A3")

load("NBselpsb.RData")
psb <- logcounts(del[hvg[1:2000],])

psb2 <- apply(psb, 1, scale)
row.names(psb2) <- colnames(psb)
psb2[1:4,1:4]
cx <- cor(t(psb2),t(psb2), method = "pearson")
cx <- cor(psb2,psb2, method = "pearson")
dim(cx)
col<- colorRampPalette(c("blue3", "white", "orangered3"))(20)

pdf("Neu3sel.pdf",  width=24, height=8, pointsize = 10)
corrplot::corrplot(cx,  col=col, order="hclust",hclust.method = "ward.D2",tl.cex=0.1)

dev.off()
pheatmap::pheatmap(cx)
