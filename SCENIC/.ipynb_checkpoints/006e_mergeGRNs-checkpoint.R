##This script merges GRNs from all the classes

suppressPackageStartupMessages({
  library(tidyverse)
})
DIR = "~/MBsnANA/HBana/SCENICana/SCENICOut/Class/"
#read metadata
setwd(DIR)
Cls <- readLines("~/MBsnANA/HBana/SCENICana/SCENIC_Scripts/RNA_Class.txt")
##This script merges GRNs from all the classes

suppressPackageStartupMessages({
  library(tidyverse)
})
DIR = "~/MBsnANA/HBana/SCENICana/SCENICOut/Class/"
#read metadata
setwd(DIR)
Cls <- readLines("~/MBsnANA/HBana/SCENICana/SCENIC_Scripts/RNA_Class.txt")
############## NEW VERSION #########################
## Based on my analysis of GRNs from old verison, I might be losing some connections
## that are indeed real and also detected in GRN of the same TF in other cell type

## 1. Obtaining a list of all GRNs-------
## First I am creating a list of GRNs from all the class GRNs. Since a TF can 
## exist in different class I am first naming each TF by adding class name in front

class_grn_list <- list()
all_tfs <- c()
for(x in Cls){
  load(paste0(x,"/SCENIC/",x,"_gmt_v2_vB.RData"))
  class_grn_list[[x]] <- regs
  tfdel <- names(regs)
  all_tfs <- rbind(all_tfs,data.frame(TF=tfdel,Class=rep(x,length(tfdel))))
  rm(regs,tfdel)
  
  
}
tfs <- unique(all_tfs$TF)
all_grns <- list()
for( x in tfs){
  del <- all_tfs[all_tfs$TF==x,"Class"]
  for(y in del){
    del_grn <- class_grn_list[[y]]
    all_grns[[paste0(x,"_",y)]] <- del_grn[[x]]
    rm(del_grn)
  }
  rm(del)
}
lengths(all_grns)
all_tfs$GRN <- paste0(all_tfs$TF,"_",all_tfs$Class)
row.names(all_tfs) <- all_tfs$GRN
save(all_grns,all_tfs, file = "AllGRNs_v2_vB.RData")
q()
#### -----
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (round(intersection/union,4))
}
overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (round(intersection/minlen,4))
}
load("AllGRNs.RData")
tfs <- names(all_grns)
row.names(all_tfs) <- paste0(all_tfs$TF,"_",all_tfs$Class)

# ##vectorized, same output as above
OS_AvB <- outer(tfs, tfs, FUN = Vectorize(function(i, j) {
  overlapS(all_grns[[i]], all_grns[[j]])
}))
rownames(OS_AvB) <- colnames(OS_AvB) <- tfs

##vectorized jaccard
JS_AvB <- outer(tfs, tfs, FUN = Vectorize(function(i, j) {
  jaccard(all_grns[[i]], all_grns[[j]])
}))
rownames(JS_AvB) <- colnames(JS_AvB) <- tfs

##ldas
colclass <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
colx <- colclass$Color
names(colx) <- row.names(colclass)
tfcolx <-sample(colorRampPalette(pal_simpsons()(14))(length(unique(all_tfs$TF))))
names(tfcolx) <- unique(all_tfs$TF)
ann_col <- list(TF=tfcolx,
                Class=colx)
ann <- all_tfs
ann$TF <- NULL
## heatmap
mb <- ((seq(1:101)-1)/100)/2
mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)
hp2 <-pheatmap(JS_AvB, 
               clustering_method = 'ward.D2',
               show_rownames = FALSE,
               show_colnames=FALSE,
               annotation_col = ann,
               annotation_colors = ann_col,
               breaks = mb, color = mc2, fontsize = 5)
pdf("TF_intersect.pdf", width = 20, height = 20, pointsize = 10)
print(hp2)
dev.off()
library(ComplexHeatmap)
library(circlize)
dmat <- as.dist(1 - JS_AvB[sel,sel]) # Using 1 - correlation as distance
hc <- hclust(dmat, method = "ward.D2") ##this look good for now
#plot(hc)
col_fun = colorRamp2(c(0,0.5, 1), c("white", "sandybrown", "brown4"))
# col_fun(seq(0,1,by=0.1))
sel_tf <- c("NFIA","PBX1","PBX3","POU2F1","RFX3","TCF12","TCF4","BACH2","NFIB")
sel <- row.names(all_tfs)[all_tfs$TF %in% sel_tf]
hp3 <- Heatmap(OS_AvB[sel,sel],
               #show_row_names = FALSE,
               show_column_names = FALSE,
               cluster_rows = hc,cluster_columns = hc,
               col=col_fun)
pdf("TF_intersect_selOS.pdf", width = 20, height = 20, pointsize = 10)
print(hp3)
dev.off()
## ----------
all_tfs_class <- list()
for(x in tfs){
  all_tfs_class[[x]] <- unique(all_tfs[all_tfs$TF==x,"Class"])
  
}
sort(lengths(all_tfs_class),decreasing = TRUE)
############### OLD VERSION ##################
## THis version focused on identifying robust GRNS by finding which TF-GRN was
## detected in both runs and then finding the robust links by finding which links
## both in terms of genes and CREs were present in both runs
## 1. Obtaining a list of all GRNs-------
## First I am creating a list of GRNs from all the class GRNs. Since a TF can 
## exist in different class I am first naming each TF by adding class name in front

class_grn_list <- list()
all_tfs <- c()
for(x in Cls){
  load(paste0(x,"/SCENIC/",x,"_gmt_v2.RData"))
  class_grn_list[[x]] <- regs
  tfdel <- names(regs)
  all_tfs <- rbind(all_tfs,data.frame(TF=tfdel,Class=rep(x,length(tfdel))))
  rm(regs,tfdel)
  
  
}
tfs <- unique(all_tfs$TF)
all_grns <- list()
for( x in tfs){
  del <- all_tfs[all_tfs$TF==x,"Class"]
  for(y in del){
    del_grn <- class_grn_list[[y]]
    all_grns[[paste0(x,"_",y)]] <- del_grn[[x]]
    rm(del_grn)
  }
  rm(del)
}
lengths(all_grns)
save(all_grns,all_tfs, file = "AllGRNs_v2.RData")
q()
#### -----
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (round(intersection/union,4))
}
overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (round(intersection/minlen,4))
}
load("AllGRNs.RData")
tfs <- names(all_grns)
row.names(all_tfs) <- paste0(all_tfs$TF,"_",all_tfs$Class)

# ##vectorized, same output as above
OS_AvB <- outer(tfs, tfs, FUN = Vectorize(function(i, j) {
  overlapS(all_grns[[i]], all_grns[[j]])
}))
rownames(OS_AvB) <- colnames(OS_AvB) <- tfs

##vectorized jaccard
JS_AvB <- outer(tfs, tfs, FUN = Vectorize(function(i, j) {
  jaccard(all_grns[[i]], all_grns[[j]])
}))
rownames(JS_AvB) <- colnames(JS_AvB) <- tfs

##ldas
colclass <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
colx <- colclass$Color
names(colx) <- row.names(colclass)
tfcolx <-sample(colorRampPalette(pal_simpsons()(14))(length(unique(all_tfs$TF))))
names(tfcolx) <- unique(all_tfs$TF)
ann_col <- list(TF=tfcolx,
                Class=colx)
ann <- all_tfs
ann$TF <- NULL
## heatmap
mb <- ((seq(1:101)-1)/100)/2
mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)
hp2 <-pheatmap(JS_AvB, 
               clustering_method = 'ward.D2',
               show_rownames = FALSE,
               show_colnames=FALSE,
               annotation_col = ann,
               annotation_colors = ann_col,
               breaks = mb, color = mc2, fontsize = 5)
pdf("TF_intersect.pdf", width = 20, height = 20, pointsize = 10)
print(hp2)
dev.off()
library(ComplexHeatmap)
library(circlize)
dmat <- as.dist(1 - JS_AvB[sel,sel]) # Using 1 - correlation as distance
hc <- hclust(dmat, method = "ward.D2") ##this look good for now
#plot(hc)
col_fun = colorRamp2(c(0,0.5, 1), c("white", "sandybrown", "brown4"))
# col_fun(seq(0,1,by=0.1))
sel_tf <- c("NFIA","PBX1","PBX3","POU2F1","RFX3","TCF12","TCF4","BACH2","NFIB")
sel <- row.names(all_tfs)[all_tfs$TF %in% sel_tf]
hp3 <- Heatmap(OS_AvB[sel,sel],
               #show_row_names = FALSE,
               show_column_names = FALSE,
               cluster_rows = hc,cluster_columns = hc,
               col=col_fun)
pdf("TF_intersect_selOS.pdf", width = 20, height = 20, pointsize = 10)
print(hp3)
dev.off()
## ----------
all_tfs_class <- list()
for(x in tfs){
  all_tfs_class[[x]] <- unique(all_tfs[all_tfs$TF==x,"Class"])
  
}
sort(lengths(all_tfs_class),decreasing = TRUE)
