##This script calculates enrichment scores (using Seurat algorithm) for a gene-set
## for hindbrain data s

suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(BiocParallel)
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
DIR_RNA <- "~/RNA"
DIR_GRN <- "~/SCENICOut/"
# 
setwd(DIR_RNA)
##1. load required data -------

# TFGRNs
load(paste0(DIR_GRN,"PN/metadata/PN_gmt.RData"))
PN_regs <- regs
names(PN_regs) <- paste0("PN:",names(PN_regs))
load(paste0(DIR_GRN,"GC/metadata/GC_gmt_v2.RData"))
GC_regs <- regs
names(GC_regs) <- paste0("GC:",names(GC_regs))
regs <- do.call(c,list(PN_regs,GC_regs))
GSEA <- regs
## 1.2 load expression data -----------
## this is compiled log transformed expression data for the entire hindbrain
load("lgcountsHB.RData")

##load metadata, combined 
load("pdPN.RData")
pdPN <- pd
load("pdGC.RData")
pdGC <- pd
pd <- rbind(pdPN,pdGC)
mdt <- mdt[,row.names(pd)]

### 1.3 remove ribosomal and mitochndiral genes from genesets--------
load("genesnoMTRibo.RData")
GSEAx <- GSEA
GSEA <- list()
cls <- names(GSEAx)
for(x in cls){
  cand <- GSEAx[[x]]
  cand <- cand[cand %in% genes]
  if(length(cand)>19){
    GSEA[[x]] <- cand
  }
  rm(cand)
}
## 2. module score----
scr <- addmodscore(mdt,features = GSEA)
colnames(scr) <- names(GSEA)
row.names(scr) <- row.names(pd)
#save(scr, pd, file = paste0("AUC_MS/pdPNGC_PNGCGRNs_MS.RData"))
save(scr, pd, file = paste0("AUC_MS/pdPNGC_RetMet_MS.RData")) ## Reactome Metabolome
#save(scr, pd, file = paste0("AUC_MS/pdPNGC_RetSig_MS.RData")) ## Reactome Metabolome
q()
## 3. AUC score----
auc <- AUCell_run(mdt, GSEA, aucMaxRank=nrow(mdt)*0.10, normAUC = FALSE,
                 BPPARAM=BiocParallel::MulticoreParam(5))
scr <- t(assay(auc)) 

save(auc, scr, pd, file =paste0("AUC_MS/pdPNGC_PNGCGRNs_AUC.RData"))
#save(auc, scr, pd, file = paste0("AUC_MS/pdPN_RetMet_AUC.RData")) ## Reactome Metabolome
save(scr, pd, file = paste0("AUC_MS/pdPNGC_RetSig_MS.RData")) ## Reactome Metabolome

q()
##### Analysis ################
setwd("~/MBsnANA/Diemana/scepro")

## 4. AUC score amrker analysis --------
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
DIR_GRN <- "~/MBsnANA/HBana/SCENICana/SCENICOut/PN/"

load("AUC_MS/pdPNGC_PNGCGRNs_MS.RData")
load("AUC_MS/pdPNGC_RetSig_MS.RData")
load("Reactome_Signaling_Pathways.RData")
rm(GSEA)
plot.data <- pd

## marker mp
cl <- unique(plot.data$Annotationlv2)
marks <- scran::findMarkers(t(scr),groups=plot.data$Annotationlv2, direction="up",
                     test.type="wilcox",pval.type="some")

sig <- c()
for(x in cl){
  del <- marks[[x]]
  del <- row.names(del)[1:3]
  sig <- rbind(sig,del)
  rm(del)
}
rm(x)
row.names(sig) <- cl
sel <- unique(c(sig))

grep("PAX",sel,value = TRUE)
write(sel, file = "metadata/TF_sel.txt")
adj <- read.delim("~/MBsnANA/HBana/SCENICana/SCENICOut/PN/metadata/")
adj <- adj[(adj$TF %in% selx) & (adj$Tar %in% selx),]
dim(adj)
write.table(adj, file ="~/MBsnANA/HBana/SCENICana/SCENICOut/PN/PN_TFadj_fil.txt",
            row.names = FALSE, quote = FALSE, sep="\t")
## 5. plot on UMAP-------
plotMS <- function(x){
  del <- zscore_c(scr[,x])
  del[del>2] <- 2
  del[del<(-2)] <- (-2)
  pdel <- cbind(plot.data,MS=del)
  pr <- ggplot(pdel %>% arrange(MS), aes(x=UMAP1, y=UMAP2 ,color=MS))+
    scale_color_gradientn(
      colours = colorRampPalette(c("deepskyblue4","grey93","magenta3"))(100),
      breaks=c(-2,0,2),
      limits=c(-2, 2))+
    ggtitle(paste0("MP",x))+
    geom_point(size=0.01)+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(pr)
}
plotMS2 <- function(x){
  del <- range_01(scr[,x])
  pdel <- cbind(plot.data,MS=del)
  pr <- ggplot(pdel %>% arrange(MS), aes(x=UMAP1, y=UMAP2 ,color=MS))+
    scale_color_gradient(low="white",high="grey4")+
    ggtitle(paste0("MP",x))+
    geom_point(size=0.01)+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(pr)
}

nx <- 8
tiff(paste0("NMF_MP",nx,"_MS.tiff"), units="in", width=4, height=4, res=300)
plotMS("BACH1")
dev.off()
## 6. modscore TF-GRN-----
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

ms2x <- apply(scr, 2, zscore_c)
#ms2x <- apply(scr, 2, range_01)

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


myColor <- colorRampPalette(c("deepskyblue4","grey93","magenta3"))(101)
# myBreaks <- c(seq(min(msd2), 0, length.out=ceiling(100/2) + 1),
#               seq(max(msd2)/100, max(msd2), length.out=floor(100/2)))
myBreaks <- c(seq(-2, 0, length.out=ceiling(100/2) + 1),
              seq(2/100, 2, length.out=floor(100/2)))
####-----
hp <- pheatmap::pheatmap(t(msd3[,sel]), 
                         clustering_method = "ward.D2",
                         #clustering_distance_cols = "correlation",
                         clustering_distance_rows = "correlation",
                         #annotation_row = cluster,
                         annotation_colors = ann_col,
                         #cluster_cols = FALSE,
                         color=myColor, breaks=myBreaks,
                         border_color = NA,
                         fontsize=5)
hp <- pheatmap::pheatmap(t(msd3[cl,sel]), 
                         clustering_method = "ward.D2",
                         clustering_distance_cols = "correlation",
                         clustering_distance_rows = "correlation",
                         cluster_cols = FALSE,
                         color=myColor, breaks=myBreaks,
                         border_color = NA,
                         fontsize = 8)
pdf("Figures/PN_Sig_MS.pdf", width = 8, height = 8,pointsize = 10)
print(hp)
dev.off()

dela <- hp$tree_row

ord_tf <- dela$labels[dela$order]

cluster = data.frame(Cl=cutree(dela, k = 3))
cluster$Cl <- paste0("Cl",cluster$Cl)
head(cluster)
colx <-sample(colorRampPalette(pal_futurama()(12))(length(unique(cluster$Cl))))
names(colx) <- unique(cluster$Cl)
ann_col <- list(Cl=colx)
write.table(cluster, file = "TFGRN_PN_hc.txt",
            sep = "\t",quote = FALSE)

delb <- dela$labels
order_g1 <- delb[dela$order]
write(order_g1, "bv2_auc.txt")
rm(dela,delb)


## 7. AUC bubble plot -----------


range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
load(paste0(DIR_GRN,"pdPN_PNGRNs_AUC.RData"))
##downsampling 500 cells
## downsampling to 1000 cells for Sig by lv2
plot.data <- c()
cl_lv <- unique(pd$Annotationlv2)
for(y in cl_lv){
  del <- pd[pd$Annotationlv2==y,]
  rn <- row.names(del)
  del2 <- del[sample(rn,1000),]
  plot.data <- rbind(plot.data,del2)
  rm(del,rn,del2)
}
rm(y)

all_cells <- row.names(plot.data)

scr2 <- scr[all_cells,]

##for reactome and signalome
colnames(scr2) <- part2[colnames(scr2),"displayName"]
scr2[1:4,1:4]

ord_mp <- mps <- part2[ord_tf,"displayName"]

cl_lv <- unique(plot.data$Annotationlv2)

auc_rank_mat <- apply(scr2, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC
table(row.names(plot.data)==row.names(auc_rank_mat))
# Set high-AUC threshold (global or per-signature)
auc_co <- 0.90

# Create result list
res_list <- list()

for (x in mps) {
  auc_vec <- auc_rank_mat[, x]
  high_vec <- auc_vec > auc_co
  
  ## create a datafrom
  df <- data.frame(cell = rownames(auc_rank_mat), 
                   high_auc = high_vec, 
                   cluster = plot.data$Annotationlv2)
  row.names(df) <- df$cell
  # Total cell above rank cut-off across clusters
  all_high <- sum(df$high_auc)
  ncells <- nrow(df)
  
  res_sig <- df %>%
    group_by(cluster) %>%
    summarise(
      tot = n(), ### total number of cells in a cluster
      a = sum(high_auc), # In cluster, high AUC
      b = tot - a,# In cluster, not high AUC
      c = all_high-a,# Outside cluster, high AUC
      d = (ncells- tot)-c,
      prop = (a / tot)*100,# % of above cut-off cells
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      odds_ratio = (a * d) / (b * c + 1e-10),
      p = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value,
      MP = x
    )
  
  res_list[[x]] <- res_sig
  rm(auc_vec,high_vec,ncells,all_high,res_sig)
}

# Combine and adjust p-values
merge_res <- bind_rows(res_list) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "BH"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-10))

merge_res <- merge_res %>%
  group_by(MP) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 5), 0)) %>%
  ungroup()

head(merge_res)
merge_res <- data.frame(merge_res)

merge_res$cluster <- factor(merge_res$cluster,levels = cl_lv)
merge_res$MP <- factor(merge_res$MP,levels = ord_mp)




p <- ggplot(merge_res, aes(x = MP, y = cluster)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 10)) +
  scale_color_gradient2(low = "white",mid = "grey87", high = "black",midpoint = 2,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        
        axis.text.x = element_text(face="bold",colour = "black",size=10,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank())+
  labs(
    x = "MetaGene programs",
    y = "Class",
    title = "Enrichment of MetaGene programs per Class",
    #subtitle = paste0("AUC > ", auc_co, " used as enrichment threshold")
  )


pdf("Figures/LDA20MP_HBlv1_MSenrprop.pdf", width = 8, height = 6, pointsize = 20)
print(p)
dev.off()
# for reactome and signlome with long names
p <- ggplot(merge_res, aes(y = MP, x = cluster)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 10)) +
  scale_color_gradient2(low = "white",mid = "grey87", high = "black",midpoint = 2,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        
        axis.text.x = element_text(face="bold",colour = "black",size=8,
                                   angle = 45,hjust=1.1,vjust = 1.1),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title=element_blank())+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) 

pdf("Figures/PNGC_topRetSigv2.pdf", width = 6, height = 10, pointsize = 10)
print(p)
dev.off()
########## Reactome ##############
## 5. marker sets ---------
load("AUC_MS/pdPN_RetSig_AUC.RData")
plot.data <- pd

## marker mp
cl <- unique(plot.data$Annotationlv3)
marks <- scran::findMarkers(t(scr),groups=plot.data$Annotationlv3, direction="up",
                            test.type="wilcox",pval.type="some")

sig <- c()
for(x in cl){
  del <- marks[[x]]
  del <- row.names(del)[1:10]
  sig <- rbind(sig,del)
  rm(del)
}
rm(x)
row.names(sig) <- cl
sel <- unique(c(sig))


## 6. modscore analysis --------

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
load("Reactome_Signaling_Pathways.RData")
load("AUC_MS/pdPN_RetSig_MS.RData")
ms2x <- apply(scr, 2, zscore_c)
msd2 <- c()

cl_lv <- unique(plot.data$Annotationlv3)
for(x in cl_lv){
  keep1 <- plot.data$Annotationlv3==x
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
colnames(msd3) <- part2[colnames(msd3),"displayName"]
selx <- part2[sel,"displayName"]
cl <- cl_lv[c(2,1,4,3,7,6,8,5)]

####-----
hp <- pheatmap::pheatmap(t(msd3[cl,selx]), 
                         #clustering_distance_cols = "correlation",
                         #clustering_distance_rows = "correlation",
                         cluster_cols = FALSE,
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
## plot on UMAP-----------
plotMS <- function(x){
  del <- zscore_c(scr[,x])
  del[del>2] <- 2
  del[del<(-2)] <- (-2)
  pdel <- cbind(plot.data,MS=del)
  pr <- ggplot(pdel %>% arrange(MS), aes(x=UMAP1, y=UMAP2 ,color=MS))+
    scale_color_gradientn(
      colours = colorRampPalette(c("deepskyblue4","grey93","magenta3"))(100),
      breaks=c(-2,0,2),
      limits=c(-2, 2))+
    ggtitle(paste0("MP",x))+
    geom_point(size=0.01)+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(pr)
}
plotMS2 <- function(x){
  del <- range_01(scr[,x])
  pdel <- cbind(plot.data,MS=del)
  pr <- ggplot(pdel %>% arrange(MS), aes(x=UMAP1, y=UMAP2 ,color=MS))+
    scale_color_gradient(low="white",high="grey4")+
    ggtitle(paste0("MP",x))+
    geom_point(size=0.01)+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "None")
  return(pr)
}

nx <- 8
tiff(paste0("NMF_MP",nx,"_MS.tiff"), units="in", width=4, height=4, res=300)
plotMS(nx)
dev.off()

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

## finding thresold using AUCell---
load("pdHBAUCCS_LDA_40comfix.RData")
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


