##This script calculates enrichment scores (using Seurat algorithm) for a gene-set
## for hindbrain data s

suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(tidyverse)
  #library(AUCell)
  #library(BiocParallel)
})
addmodscore <- function(obj, features, name = "Cluster",
                        pool = row.names(obj), nbin=24,
                        ctrl = 100, k =FALSE, seed = 123){
  # Find how many gene lists were provided. In this case just one.
  cluster.length <- length(x = features)
  # For all genes, get the average expression across all cells (named vector)
  data.avg <- Matrix::rowMeans(x = obj[pool, ])
  # Order genes from lowest average expression to highest average expression
  data.avg <- data.avg[order(data.avg)]
  
  # Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations.
  # The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
  data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                  n = nbin,
                                  labels = FALSE,
                                  right = FALSE)
  # Set the names of the cuts as the gene names
  names(x = data.cut) <- names(x = data.avg)
  # Create an empty list the same length as the number of input gene sets.
  # This will contain the names of the control genes
  ctrl.use <- vector(mode = "list", length = cluster.length)
  
  # For each of the input gene lists:
  for (i in 1:cluster.length) {
    # Get the gene names from the input gene set as a character vector
    features.use <- features[[i]]
    
    # Loop through the provided genes (1:num_genes) and for each gene,
    # find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
    for (j in 1:length(x = features.use)) {
      # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number.
      # We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                          size = ctrl,
                                          replace = FALSE)))
    }
  }
  # Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin,
  # there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling
  # the same gene more than once.
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  
  
  ## Get control gene scores
  
  # Create an empty matrix with dimensions;
  # number of rows equal to the number of gene sets (just one here)
  # number of columns equal to number of cells in input Seurat obj
  ctrl.scores <- matrix(data = numeric(length = 1L),
                        nrow = length(x = ctrl.use),
                        ncol = ncol(x = obj))
  
  # Loop through each provided gene set and add to the empty matrix the mean expression of the control genes in each cell
  for (i in 1:length(ctrl.use)) {
    # Get control gene names as a vector
    features.use <- ctrl.use[[i]]
    # For each cell, calculate the mean expression of *all* of the control genes
    ctrl.scores[i, ] <- Matrix::colMeans(x = obj[features.use,])
  }
  
  ## Get scores for input gene sets
  # Similar to the above, create an empty matrix
  features.scores <- matrix(data = numeric(length = 1L),
                            nrow = cluster.length,
                            ncol = ncol(x = obj))
  
  # Loop through input gene sets and calculate the mean expression of these genes for each cell
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- obj[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  
  # Subtract the control scores from the feature scores -
  # the idea is that if there is no enrichment of the genes in the geneset in a cell,
  # then the result of this subtraction should be ~ 0
  features.scores.use <- features.scores - ctrl.scores
  
  # Name the result the "name" variable + whatever the position the geneset was in the input list, e.g. "Cluster1"
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  
  # Change the matrix from wide to long
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  
  # Give the rows of the matrix, the names of the cells
  rownames(x = features.scores.use) <- colnames(x = obj)
  
  return(features.scores.use)
  
}

##Directories 
Cls <- commandArgs(trailingOnly = TRUE)
DIR_RNA <- "~/MBsnANA/Diemana/scepro/"
setwd(DIR_RNA)

##1. load required data -------
##TF-GRNs
load("~/MBsnANA/HBana/SCENICana/SCENICOut/Class/AllGRNs_v2.RData")
grns <- row.names(all_tfs) <- paste0(all_tfs$TF,"_",all_tfs$Class)
sel_grns <- grns[all_tfs$Class==Cls]
#write(grns, file="~/MBsnANA/HBana/SCENICana/SCENICOut/Class/AllGRNsnames.txt")

GSEA <- list()
#GSEA[[grn]] <-  all_grns[[grn]]
for(x in sel_grns){
  GSEA[[x]] <- all_grns[[x]]
}

## this is compiled log transformed expression data for the entire hindbrain
load("lgcountsHB.RData")

##load metadata, combined
load("plotdataHBTempANN250305.RData")
plot.data <- plot.data[plot.data$Annotationlv3!="ND",]

#for fast run subset plotdata
# cells <- row.names(plot.data)
# cells_sel <- sample(cells,0.25*length(cells))
# plot.data <- plot.data[cells_sel,]

mdt <- mdt[,row.names(plot.data)]
# # ##module score----
scr <- addmodscore(mdt,features = GSEA)
colnames(scr) <- names(GSEA)
row.names(scr) <- row.names(plot.data)

save(scr, plot.data, file = paste0("AUC_MS/pdHB_",Cls,"_MS.RData"))

rm(scr)
q()
## merge scores from all splits ------
setwd("~/MBsnANA/Diemana/scepro/")
zscore_c <- function(x){
  sdx <- sd(x)
  mx <- mean(x)
  ot <- boxplot.stats(x)$out
  q <- boxplot.stats(x)$stats
  keep <- (x < q[1])
  x[keep] <- q[1]
  keep <- (x > q[5])
  x[keep] <- q[5]
  x2 <- x-mx
  z <- x2/sdx
  return(z)
}

Cls <- readLines("~/MBsnANA/Diemana/HB_Scripts/RNA_Class.txt")

scr_com <- c()
scr_com_z <- c()

for(x in Cls){
  load(paste0("AUC_MS/pdHB_",x,"_MS.RData"))
  #scrx <- apply(scr, 2, zscore_c)
  scr_com <- cbind(scr_com,as.matrix(scr))
  #scr_com_z <- cbind(scr_com_z,scrx)
  rm(scr)
}
rm(x)
## i added scr_com (un-scaled) and renamed scr_com to scr_com_z on 180825
save(scr_com,scr_com_z, plot.data, file = "AUC_MS/pdHB_AllGRNs_MSz.RData") 

## ------------
setwd("~/MBsnANA/Diemana/scepro/")
load("AUC_MS/pdHB_AllGRNs_MSz.RData")

load("pdHB_LDAMP_40comfix_100g_AUC.RData")
mps <- names(cells_assignment)

thr <- data.frame(MPs=names(cells_assignment))
row.names(thr) <- thr$MPs
thrs <- c()
for(x in mps){
  del <- cells_assignment[[x]]
  del<- del$aucThr$thresholds
  del_dim <- dim(del)[1]
  if(del_dim==3){
    del <- rbind(del,c(0,0))
    row.names(del)[4] <- "minimumDens"
  }
  thrs <- rbind(thrs,del[,1])
  rm(del,del_dim)
}
rm(x)
row.names(thrs) <- mps
thrs[1:4,]
thrs <- data.frame(thrs)
thrs$Selected <- 0

scr[1:4,1:4]
plist <- list()
for(x in mps){
  del_df <- data.frame(Score=scr[,x])
  a <-thrs[x,"Global_k1"]
  b <-thrs[x,"L_k2"]
  c <-thrs[x,"R_k3"]
  d <-thrs[x,"minimumDens"]
  e <-thrs[x,"Selected"]
  pdel <- ggplot(del_df,aes(x=Score))+
    geom_histogram()+
    scale_fill_manual(values = c("grey78"))+
    theme_bw()+
    geom_vline(xintercept = c(a,b,c,d,e),
               color=c("red","blue","green","orange","brown"),
               linetype=c(rep("dashed",4),"solid"))+
    ggtitle(x)
  plist[[x]] <- pdel
  rm(pdel,del_df)

}
rm(x)
pdf("LDA_AUCthr_40comfix.pdf",width = 20,height = 10,pointsize = 10)
gridExtra::grid.arrange(grobs=plist,ncol=8,nrow=5)
dev.off()
write.table(thrs, file = "LDA_AUCthresold.txt",
            sep="\t",quote = FALSE)
thrs <- read.delim("LDA_AUCthresold.txt",row.names = 1)
####  range 0 to 1-----
ms2x <- apply(scr, 2, range_01)
msd2 <- c()

cl_lv <- unique(plot.data$Annotationlv2)
for(x in cl_lv){
  keep1 <- plot.data$Annotationlv2==x
  msx <- ms2x[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]
msd3 <- apply(msd2, 2, range_01)
msd3[1:4,1:4]

mc <- colorRampPalette(c("white","grey75","black"))(101)
mb <- seq(1:101)-1
mb <- mb/100
hp <- pheatmap::pheatmap(msd3, 
                         clustering_method = "ward.D2",
                         clustering_distance_cols = "correlation",
                         clustering_distance_rows = "correlation",
                         #cluster_rows = FALSE,
                         color=mc, breaks=mb,
                         border_color = NA,fontsize = 8)

pdf("Figures/NMFMP40HBannlv2enr01_nc.pdf", width = 8, height = 16, pointsize = 10)
print(hp)
dev.off()

#myColor <- colorRampPalette(c(rev(brewer.pal(9,"PuBu")),"white",brewer.pal(9,"OrRd")))(100)
#myBreaks <- ((seq(1:101)-1)/100)
####  zscore -----
ms2x <- apply(scr, 2, zscore_c)
msd2 <- c()

cl_lv <- unique(plot.data$Annotationlv2)
for(x in cl_lv){
  keep1 <- plot.data$Annotationlv2==x
  msx <- ms2x[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]
msd3 <- apply(msd2,2,zscore_c)
msd3[1:4,1:4]

mc <- colorRampPalette(c("deepskyblue4","grey93","magenta3"))(101)
# mb <- c(seq(min(msd2), 0, length.out=ceiling(100/2) + 1),
#               seq(max(msd2)/100, max(msd2), length.out=floor(100/2)))
mb <- c(seq(-2, 0, length.out=ceiling(100/2) + 1),
              seq(2/100, 2, length.out=floor(100/2)))
####-----
hp <- pheatmap::pheatmap(msd2, 
                         clustering_method = "complete",
                         #clustering_distance_cols = "correlation",
                         #clustering_distance_rows = "correlation",
                         #cluster_rows = FALSE,
                         color=mc, breaks=mb,
                         border_color = NA,
                         fontsize = 8)
pdf("Figures/MarkersHBannlv2enrZAUC.pdf", width = 16, height = 16, pointsize = 10)
print(hp)
dev.off()

dela <- hp$tree_row
tiff( "tree_bv2_AUC.tiff", width = 18, height = 6, units = "in", res = 300)
plot(dela, cex=0.5)+abline(h=c(0.4,0.5,0.8), col="red")
dev.off()

delb <- dela$labels
order_g1 <- delb[dela$order]
write(order_g1, "bv2_auc.txt")
rm(dela,delb)
### enrichment by 0 to 1 scaling -------
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
load("pdHB_NMFMP20_AUC.RData")
table(row.names(plot.data) %in% row.names(scr))
scr <- scr[row.names(plot.data),]
mdt2 <- apply(scr,2,range_01)
mdt2[1:4,1:4]
annlv2 <- unique(plot.data$Annotationlv2)


msd <- c()
for(x in annlv2){
  cells <- row.names(plot.data)[plot.data$Annotationlv2==x]
  scrx <- colMeans(mdt2[cells,])
  msd <- rbind(msd,scrx)
  rm(cells,scrx)
  
}
row.names(msd) <- annlv2
msd[1:4,1:4]



