## plot seqlet enrichement in marker CREs using homer based analysis
## analysizng module or AUC score for MAPS and SASP in normal and tumors
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggrepel)
})

setwd("~/ATACana/Outs/")

## finding enrichment of sqlets in class/cluster based on FDR and proportion

cls <- readLines("~/RNA/RNA_Class.txt")
##these are syntaxes/motifs identified from DeepHB analysis
## each motif is assigned an id randing from seqid_001 to seqid_165
seqids <- readLines("~/DeepHB/Outs/seqids.txt")
##these are selected motifs that are plotted
sel_seq <-  paste0("seqid_",sprintf("%03d", c(5,1,3,17,20,133))) ##for class level

## read pvalue and proportion data from HOMER
DIR <- "~/ATACana/Outs/markerCREs/motifenr/known_tab/"
seq_prop <- c()
seq_pval <- c()
for(x in cls){
    del <- read.delim(paste0(DIR,x,".txt"),header=FALSE,row.names=1)
    del <- del[seqids,]
    seq_pval <- cbind(seq_pval,del$V3)
    del$V7 <- as.numeric(gsub("%", "", del$V7))
    seq_prop <- cbind(seq_prop,del$V7)
    rm(del)
    }
colnames(seq_pval) <- colnames(seq_prop) <-cls
row.names(seq_pval) <- row.names(seq_prop) <-seqids
seq_pval[1:4,1:4]
seq_prop[1:4,1:4]


##fdr correct pvalues
seq_fdr <- matrix(
  p.adjust(as.vector(seq_pval), method = "fdr"),
  nrow = nrow(seq_pval),
  ncol = ncol(seq_pval)
)
colnames(seq_fdr) <- colnames(seq_pval)
row.names(seq_fdr) <- row.names(seq_pval) 
seq_fdr <- -log10(seq_fdr)
seq_fdr[seq_fdr<2] <-0
seq_fdr <- round(seq_fdr,2)
seq_fdr[1:4,1:4]

## merge FDR and prop data
dd1 <- reshape2::melt(seq_fdr[sel_seq,])
dd2 <- reshape2::melt(seq_prop[sel_seq,])
row.names(dd1) <- paste0(dd1$Var1,"_",dd1$Var2)
row.names(dd2) <- paste0(dd2$Var1,"_",dd2$Var2)
table(row.names(dd1)==row.names(dd2))
head(dd1)
head(dd2)
colnames(dd1) <- c("Seqid","Class","FDR")
dd1$Prop <- dd2$value
dd1$FDR[dd1$FDR>100] <- 100

dd1$Class <- factor(dd1$Class,levels=cls)
dd1$Seqid <- factor(dd1$Seqid,levels=sel_seq)


dd1x <- dd1[dd1$Seqid %in% sel_seq,]
dd1x$Seqid <- factor(dd1x$Seqid,levels=sel_seq)
p <- ggplot(dd1x, aes(y = Seqid, x= Class)) +
  geom_point(aes(size = FDR,color = Prop)) +
scale_color_gradient2(low = "white",mid = "grey70", high = "black",midpoint = max(dd1x$Prop)/2,
                        name = "Proportion")+
  scale_size_continuous(range = c(0, 10)) +
    theme_light()+
  theme(axis.line = element_blank(),
        panel.border=element_rect(colour = 'black',linewidth=1),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.x = element_text(colour = "black",size=10,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(face="italic",colour = "black",size=10),
        axis.title=element_blank())+
  scale_y_discrete(position = "right")

## for progenitor
cls <- readLines("Pro_clusters.txt")
seqids <- readLines("~/DeepHB/Outs/seqids.txt")
DIR <- "~/ATACana/Outs/markerCREs/Pro_cluster/motifenr/known_tab/"
seq_prop <- c()
seq_pval <- c()
for(x in cls){
    del <- read.delim(paste0(DIR,x,".txt"),header=FALSE,row.names=1)
    del <- del[seqids,]
    seq_pval <- cbind(seq_pval,del$V3)
    del$V7 <- as.numeric(gsub("%", "", del$V7))
    seq_prop <- cbind(seq_prop,del$V7)
    rm(del)
    }
colnames(seq_pval) <- colnames(seq_prop) <-cls
row.names(seq_pval) <- row.names(seq_prop) <-seqids
seq_pval[1:4,1:4]
seq_prop[1:4,1:4]
dd1$Class <- factor(dd1$Class,levels=rev(cls))
###these ids are plottted
sel_seq <-  paste0("seqid_",sprintf("%03d", c(36,66,12,34,5,20,26)))
dd1x <- dd1[dd1$Seqid %in% sel_seq,]
dd1x$Seqid <- factor(dd1x$Seqid,levels=sel_seq)
p <- ggplot(dd1x, aes(y = Seqid, x= Class)) +
  geom_point(aes(size = FDR,color = Prop)) +
scale_color_gradient2(low = "white",mid = "grey70", high = "black",midpoint = max(dd1x$Prop)/2,
                        name = "Proportion")+
  scale_size_continuous(range = c(0, 10)) +
    theme_light()+
  theme(axis.line = element_blank(),
        panel.border=element_rect(colour = 'black',linewidth=1),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.x = element_text(colour = "black",size=14,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title=element_blank())




