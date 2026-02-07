## This script extract Pontine nuclei lineage from hindbrain
suppressPackageStartupMessages({
  library(Matrix)
  library(ggsci)
  library(ggrepel)
  library(slingshot)
  library(tidyverse)
  library(scran)
})
DIR <- "~/RNA/"
setwd(DIR)
## load data---------
## load plotdata
load("plotdataHBv2.RData")
cl <- c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES",
        "Neurons_I:PonNuc_def","Neurons_I:PonNuc_diff_early","Neurons_I:PonNuc_diff_late",
        "Neurons_II:PonNuc_mat")
pd <- plot.data[plot.data$Cluster %in% cl,]
rm(plot.data)

## load NMF
load("NMF_HB.RData") ##this is NMF factors from entire hindbrain
load("plotdataHB") ##fixed with old and CB annotation
dim(rd_nmf)
row.names(rd_nmf) <- row.names(plot.data)
colnames(rd_nmf) <- paste0("NMF_",seq(100))
nmf <- rd_nmf[row.names(pd),]
pd$Stage <- plot.data[row.names(pd),"Stage"]
rm(plot.data,rd_nmf)

## constrcut UMAP
set.seed(145)
del <- uwot::umap(nmf,min_dist = 0.3,metric = "cosine")
pd$UMAP1 <- del[,1]
pd$UMAP2 <- del[,2]
## get HVG HEG an dTf -----------
save(pd, file = "pdPN.RData")

## obtian logcounts
load("lgcountsHB.RData")
gsm <- readLines("~/extrafiles/genesGSM.txt")
pc <- read.table("~/extrafiles/gencdH38p13r37_PC.txt", header = TRUE)
pc <- as.character(pc$Gene.ID)
gsm <- gsm[gsm %in% pc]
gsm <- gsm[gsm %in% row.names(mdt)]
mdt <- mdt[gsm,row.names(pd)]

## finding TFs that are expressed in at least 10% of the cells in a cluster
tfs <- readLines("~/SCENICfiles/allTFs_hg38.txt")
tfs <- tfs[tfs %in% row.names(mdt)]
mdt_tf <- mdt[tfs,]

exp <- c()
annlv2 <- unique(pd$Annotationlv3)
for(x in annlv2){
  mdt_tfx <- mdt_tf[,pd$Annotationlv3==x]
  del_mat <- as(matrix(0,dim(mdt_tfx)[1],dim(mdt_tfx)[2]),"dgCMatrix")
  del_mat[mdt_tfx>0]<-1
  del <- rowSums(del_mat)
  del <- (del/dim(mdt_tfx)[2])*100
  exp <- cbind(exp,del)
  rm(mdt_tfx,del_mat,del)
}
rm(x)
row.names(exp) <- tfs
colnames(exp) <- annlv2

#exp[1:4,1:2]
keep <- apply(exp,1,function(x){any(x >=50)})
table(keep)

tfs_sel <- tfs[keep]
##downsampling large annotations to max 3k cells
tb <- data.frame(table(pd$Annotationlv3))
cl <- as.character(tb$Var1)[tb$Freq>1500]
keep <- pd$Annotationlv3 %in% cl
table(keep)
pd1 <- pd[keep, ]
pd2 <- pd[!keep, ]
pd3 <- c()
rm(keep,tb)
for(y in cl){
  del <- pd[pd$Annotationlv3==y,]
  rn <- row.names(del)
  del2 <- del[sample(rn,1500),]
  pd3 <- rbind(pd3,del2)
  rm(del,rn,del2)
}
rm(y)
pd <- rbind(pd2,pd3)
table(pd$Annotationlv2)
table(pd$Annotationlv3)
save(pd, file = "pdPNdwn.RData")

# find HEG
gene_means <- rowMeans(mdt)
gene_means <- sort(gene_means, decreasing = TRUE)
heg <- names(gene_means)
load("~/extrafiles/genesnoRPMT.RData")
heg <- heg[heg %in% genes]

#not removing sex genes
# load("~/MBATAC/MBnewana/HGNCsexchr.RData")
# heg <- heg[!(heg %in% SX)]
heg <- heg[1:500]

## finding HVGs
dec <- modelGeneVar(mdt)
dec <- data.frame(dec)
dec <- dec[!is.nan(dec$FDR),]
dec<- dec[dec$bio>0,]
dec <- dec[dec$mean>0.02,]
dec <- dec[order(dec$bio,decreasing = TRUE),]
hvg <- row.names(dec)
hvg <- hvg[hvg %in% genes]

save(hvg,heg,tfs_sel, file = "~/SCENICOut/PN/PN_HVG_pyscenic.RData")
q()
