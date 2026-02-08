## This script performs analysis of context-specific TF-GRN activity across
## classes. 
## First I plot a correlation heatmap, then perform pairwise wilcoxon to identify
## signficance of both TF-GRN gene activity score and associated CRE activity score
## 1. correlation heatmap -------
## For TF associated with multiple classes generate correlation values and order of classes

suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(pheatmap)
  library(scran)
  library(scater)
  
})
setwd("~/RNA/")
tf <- commandArgs(trailingOnly = TRUE) ##selected TF
## 1.a GRN sets -------
load("~/SCENIC/AllGRNs.RData")

all_tfs <- all_tfs[all_tfs$TF==tf,]
sel_cls <- all_tfs$Class
sel_grns <- all_tfs$GRN
## defined order for two selected TFs
if(tf=="NFIC"){
  sel_cls <- c("GC","Neurons_I","Neuroblasts_III","CBINT",
               "MixedGlia","MES_VAS")
  sel_grns <- paste0(tf,"_",sel_cls)
}else if(tf=="ARNT2"){
  sel_cls <- c("Neurons_I","Neurons_II","Neuroblasts_II","Neuroblasts_III",
               "PreO","Progenitors","MixedGlia","Glioblasts")
  sel_grns <- paste0(tf,"_",sel_cls)
}

## 1.b module score of each TF-GRNs ----------
load("AUC_MS/pdHB_AllGRNs_MSz.RData")
kk <- colnames(scr_com)

load("plotdataHBsv2.RData")
table(row.names(plot.data)==row.names(scr_com_z))
plot.data <- plot.data[plot.data$Class %in% sel_cls,]


scr_com_z <- scr_com_z[row.names(plot.data),sel_grns]
scr_com <- scr_com[row.names(plot.data),sel_grns]

## 1.c correlation -----
corx <- cor(scr_com_z,scr_com_z)
hp <- pheatmap(corx,clustering_method = "ward.D2",
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               silent = TRUE)

ord <- hp$tree_row
ord <- ord$labels[ord$order]
## fix some order, this was already defined above
if((tf %in% c("ARNT2","NFIC"))){
  ord <- sel_grns
} else {
  sel_grns <- ord
}
ord2 <- all_tfs[ord,"Class"]
write(ord2, file = paste0("TFcontext/",tf,"_ordofclass.txt"))

corx <- corx[ord,ord]
colnames(corx) <- row.names(corx) <- ord2
write.table(corx, file = paste0("TFcontext/",tf,"_grn_cor.txt"),sep = "\t",quote = FALSE)
rm(scr_com_z)


## 2. RNA-part: plot TF-GRN AUC signficance test -----------
## for this subset data
ord2 <- readLines(paste0("TFcontext/",tf,"_ordofclass.txt"))
ord <- paste0(tf,"_",ord2)
## for TF detected in more than 7 classes, the memory requirements for ATAC step is
## quite high, and hence i subsetted ATAC data.
## to be consistent I also RNA data to  max 6000 cells per class
cellsx <- readLines("CellHBRNAsubsetclass.txt")
cells <- row.names(plot.data)
cells <- intersect(cells, cellsx)
plot.data <- plot.data[cells,]
scr_com <- scr_com[cells,]
## 2.1 signficance, making bubble plot based on AUC and FDR

## using pair-wise wilcox and finding average values
## The idea is to first perform pairwise comparison of each class to other classes
## one by one for each of the gene-set. Then I generate a combined data
marks <- pairwiseWilcox(t(scr_com), groups=plot.data$Class,direction="up")
sets <- data.frame(A=marks$pairs$first,B=marks$pairs$second)
dd <- marks$statistics
names(dd) <- paste0(sets$A," vs ",sets$B)

pval_df <- c()

for(x in ord2){
  b <- sets[sets$A==x,"B"]
  df_del <- c()
  for(y in b){
    z <- paste0(x," vs ",y)
    del <- data.frame(dd[[z]])
    del$B <-y
    del$GRN <- row.names(del)
    row.names(del) <- NULL
    df_del <- rbind(df_del,del[,c("GRN","B","AUC","p.value")])
    rm(del,z)
  }
  rm(y)
  
  #before finding fisher method adjusted p-value, fix the 0 pvalue
  ## either to the next minimum value or 1e-100, which ever is smaller
  min_p <- min(df_del$p.value[df_del$p.value!=0])
  if(min_p<1e-100){
    df_del[df_del$p.value==0,"p.value"] <- min_p
  }else{
    df_del[df_del$p.value==0,"p.value"] <- 1e-100
  }
  df_del$Class <- x
  pval_df <- rbind(pval_df,df_del)
  rm(min_p,df_del,b)
}
rm(x)
head(pval_df) ## this is the combined data

## now from the combined data, I want to identify for each GRN which classes are significant
## so I have to extract GRNs and each class. this will have data for n-1 comparisons.
## from that  I identify mean AUC and a combined pvalue (fisher method)
## then I want to see if that GRN is signficantly enriched in that class by checking
## if in majority of comparisons, that class is enriched
## i am checking signficance AUC>0.51 and fdr<1e-5, then I check the proprtion
## of comparisons the class is enriched has to be more than >=0.5

## 2.2 merge p-values across multiple comparisons

## this is to combine pvalue from various pair-wise test
combine_fisher <- function(p) {
  stat <- -2 * sum(log(p))
  stats::pchisq(stat, df = 2*length(p), lower.tail = FALSE)
}

alpha <- 1e-5
pd <- pval_df %>%
  dplyr::group_by(GRN, Class) %>%
  dplyr::summarise(
    mean_AUC = mean(AUC, na.rm = TRUE),## mean AUC of a GRN in a class
    prop_sig = mean((AUC > 0.51) & (p.adjust(p.value, method = "bonferroni") < alpha)),## proportion of significant enrcihments
    fisher_p   = combine_fisher(p.value),##combined pvalue for this GRN for this Class
    Sig = ifelse(prop_sig > .5,"Yes","No"),## is it signficantly enriched?
    .groups = "drop"
  )

pd <- data.frame(pd)
head(pd)
pd[pd$mean_AUC<.1,"mean_AUC"] <- .1
pd$FDR <- p.adjust(pd$fisher_p,method = "bonferroni") # final correction of fisher_p
pd[pd$FDR==0,"FDR"] <- min(pd$FDR[pd$FDR!=0])/10 ## take care of zeros
#pd[pd$FDR<1e-10,"FDR"] <- 1e-10
pd$FDR <- -log10(pd$FDR)
pd[pd$mean_AUC<0.5,"Sig"] <- "No" ## if mean_AUC was <0.5 then remove from signficance

write.table(pd, file = paste0("TFcontext/",tf,"_pvaldf_MS.txt"),sep = "\t",quote = FALSE,row.names = FALSE)

rm(scr_com,plot.data,pd)

## 3. ATAC-part: plot TF-GRN AUC signficance test-----------
suppressPackageStartupMessages({
  library(ArchR)
  library(AUCell)
})
DIR_GRN <- "~/SCENIC/Class/"
DIR_ATAC <- "~/ATACana/Outs/"

## 3.1 extract region sets
## Only CREs associated with a TF in a class are extracted
regCREs <- c()


for(i in 1:length(ord2) ){
  x <- ord2[i]
    y <- ordx[i]
  df <- read.delim(paste0(DIR_GRN,y,"/SCENIC/",y,"_eRegulons_fix.txt"))
  df <- df[df$TF==tf,c("TF","Region","Gene")]
  df$GRN <- ord[i]
  df$Class <- x
  regCREs<- rbind(regCREs,df)
  rm(df)
  
}
rm(i)
head(regCREs)
regCREs$X <- paste0(regCREs$Class,"_",regCREs$Region)
## duplciate CREs for chosen TF in a particular class are removed
regCREs <- regCREs[!duplicated(regCREs$X),]
row.names(regCREs) <- regCREs$X
peaks <- unique(regCREs$Region)

## cre-set for AUC
GSEA <- list()
for(x in ord2){
  GSEA[[x]] <- regCREs[regCREs$Class==x,"Region"]
}
lengths(GSEA)


## 3.2 load ATAC metadata
load("~/ATACana/Outs/plotdataATACv6.RData")
## for TF detected in more than 7 classes, the memory requirements are
## quite high to perform below analysis, so I subsetted it to max 4000 cells per class
##subset 
cell_atac <- readLines("~/ATACana/Outs/CellHBATACsubsetclass.txt")
plot.data <- plot.data[cell_atac,]
plot.data <- plot.data[plot.data$Class%in%ord2,]
cells <- row.names(plot.data)
dim(plot.data)

## 3.2 load Peak matrix for AUC score
proj <- loadArchRProject(path = paste0(DIR_ATAC,"comATACx/"))

##extract PeakMatrix for selected cells
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)

## get Matrix and subset to peaks
region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)

PM <- assay(PeakMatrix)
row.names(PM) <- row.names(region_names)
PM[1:4,1:4]

## since calculating AUC could be a huge memeory issue for entire
## dataset, I am doing it per class and then concatenating
scr <- c()
for(x in ord2){
  cells_sel <- row.names(plot.data)[plot.data$Class==x]
  PMx <- PM[,cells_sel] ###ver important to make sure cell order is same
  ms <- AUCell_run(PMx, GSEA, aucMaxRank=nrow(PM)*0.10, normAUC = FALSE,
                   BPPARAM=BiocParallel::MulticoreParam(6))
  scr_del <- t(assay(ms))
  colnames(scr_del) <- paste0(tf,"_",colnames(scr_del))
  scr <- rbind(scr,scr_del)
  rm(cells_sel,PMx,ms,scr_del)
}
scr <- scr[row.names(plot.data),]


#save(plot.data,scr, file =paste0("AUC_MS/pdAUC_",tf,"_CRE.RData")) #TF GRN CREs

## 3.4 signficance, making bubble plot based on AUC and FDR

## using pair-wise wilcox and finding average values
marks <- pairwiseWilcox(t(scr), groups=plot.data$Class,direction="up")
sets <- data.frame(A=marks$pairs$first,B=marks$pairs$second)
dd <- marks$statistics
names(dd) <- paste0(sets$A," vs ",sets$B)

pval_df <- c()

for(x in ord2){
  b <- sets[sets$A==x,"B"]
  df_del <- c()
  for(y in b){
    z <- paste0(x," vs ",y)
    del <- data.frame(dd[[z]])
    del$B <-y
    del$GRN <- row.names(del)
    row.names(del) <- NULL
    df_del <- rbind(df_del,del[,c("GRN","B","AUC","p.value")])
    rm(del,z)
  }
  rm(y)
  
  #before finding fisher method adjusted p-value, fix the 0 pvalue
  ## either to the next minimum value or 1e-100, which ever is smaller
  min_p <- min(df_del$p.value[df_del$p.value!=0])
  if(min_p<1e-100){
    df_del[df_del$p.value==0,"p.value"] <- min_p
  }else{
    df_del[df_del$p.value==0,"p.value"] <- 1e-100
  }
  df_del$Class <- x
  pval_df <- rbind(pval_df,df_del)
  rm(min_p,df_del,b)
}
rm(x)
head(pval_df)

## this is to combine pvalue from various pair-wise test
combine_fisher <- function(p) {
  stat <- -2 * sum(log(p))
  pchisq(stat, df = 2*length(p), lower.tail = FALSE)
}

alpha <- 1e-5
pd <- pval_df %>%
  dplyr::group_by(GRN, Class) %>%
  dplyr::summarise(
    mean_AUC = mean(AUC, na.rm = TRUE),## mean AUC of a GRN in a class
    prop_sig = mean((AUC > 0.51) & (p.adjust(p.value, method = "bonferroni") < alpha)),## proportion of significant enrcihments
    fisher_p   = combine_fisher(p.value),##combined pvalue for this GRN for this Class
    Sig = ifelse(prop_sig > .5,"Yes","No"),## is it signficantly enriched?
    .groups = "drop"
  )

pd <- data.frame(pd)
head(pd)
pd[pd$mean_AUC<.1,"mean_AUC"] <- .1

pd$FDR <- p.adjust(pd$fisher_p,method = "bonferroni") # final correction of fisher_p
pd[pd$FDR==0,"FDR"] <- min(pd$FDR[pd$FDR!=0])/10 ## take care of zeros
#pd[pd$FDR<1e-10,"FDR"] <- 1e-10
pd$FDR <- -log10(pd$FDR)
pd[pd$mean_AUC<0.5,"Sig"] <- "No" ## if mean_AUC was <0.5 then remove from signficance


write.table(pd, file = paste0("TFcontext/",tf,"_pvaldf_ATAC.txt"),sep = "\t",quote = FALSE,row.names = FALSE)


q()

## 4. Plots maps ---------------
# ## example plot
ord2 <- readLines(paste0("TFcontext/",tf,"_ordofclass.txt"))
corx <- read.delim(file=paste0("TFcontext/",tf,"_grn_cor.txt"),header = TRUE,row.names = 1)
corx$Class <- row.names(corx)
pd_corx <- reshape2::melt(corx)
pd_corx$size <- abs(pd_corx$value)

pd_corx$Class <- factor(pd_corx$Class, levels = ord2)
pd_corx$variable <- factor(pd_corx$variable, levels = rev(ord2))

p <- ggplot(pd_corx, aes(x = Class, y = variable, fill = value,size = size)) +
  geom_point(shape = 21) +
  scale_size(range = c(0, 10),limits = c(0,1),  name = "Absolute correlation")+
  scale_fill_gradient2(
    low = "blue3", mid = "white", high = "orangered2",
    midpoint = 0, limits = c(-1, 1),name = "Correlation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color="black",size=14),
    panel.grid.major = element_line(color="grey50",linewidth = .5),
    panel.grid.minor = element_blank()
  )+
  scale_x_discrete(position = "top")




## plot AUC maps for a TF's GRN across all the classes ----------
pd <- read.delim(paste0("TFcontext/",tf,"_pvaldf_MS.txt"),header = TRUE)
pd$Sig <- factor(pd$Sig,levels=c("No","Yes"))
pd$GRN <- paste0("g",pd$GRN)
pd$GRN <- factor(pd$GRN, levels = rev(paste0("g",tf,"_",ord2)))

pd$Class <- factor(pd$Class, levels = ord2)


p<- ggplot(pd, aes(y = GRN, x = Class,color=Sig,
                   size = mean_AUC, fill = mean_AUC)) +
  geom_point(shape = 21) +
  scale_color_manual(values = c("grey87","black"))+
  #scale_size(range = c(2, 12),  name = "Log10FDR") +
  scale_size(range = c(1, 12),limits = c(0.1,1),  name = "Mean AUC") +
  #scale_fill_viridis_c(option = "blues",direction = -1, name = "Mean AUC") +
  scale_fill_gradientn(colors = c("white", "lightblue", "turquoise4"),
                       limits = c(0.1, 1),
                       name = "Mean AUC") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(face="italic"),
    axis.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color="black",size=14),
    panel.grid.major = element_line(color="grey50",linewidth = .5),
    panel.grid.minor = element_blank()
  )+
  scale_x_discrete(position = "top")

w=.27*length(ord2)+4
h=.2*length(ord2)+3


## plot AUC maps CRE for a TF's GRN across all the classes ----------

pd_ATAC <- read.delim(paste0("TFcontext/",tf,"_pvaldf_MS.txt"),header = TRUE)
pd_ATAC$Sig <- factor(pd_ATAC$Sig,levels=c("No","Yes"))
pd_ATAC$GRN <- paste0("g",pd_ATAC$GRN)
pd_ATAC$GRN <- factor(pd_ATAC$GRN, levels = rev(paste0("g",tf,"_",ord2)))

pd_ATAC$Class <- factor(pd_ATAC$Class, levels = ord2)

p<- ggplot(pd, aes(y = GRN, x = Class,color=Sig,
                   size = mean_AUC, fill = mean_AUC)) +
  geom_point(shape = 21) +
  scale_color_manual(values = c("grey87","black"))+
  #scale_size(range = c(2, 12),  name = "Log10FDR") +
  scale_size(range = c(1, 12),limits = c(0.1,1),  name = "Mean AUC") +
  #scale_fill_viridis_c(option = "blues",direction = -1, name = "Mean AUC") +
  scale_fill_gradientn(colors = c("white", "lightblue", "turquoise4"),
                       limits = c(0.1, 1),
                       name = "Mean AUC") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(face="italic"),
    axis.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color="black",size=14),
    panel.grid.major = element_line(color="grey50",linewidth = .5),
    panel.grid.minor = element_blank()
  )+
  scale_x_discrete(position = "top")

w=.27*length(ord2)+4
h=.2*length(ord2)+3

q()
