## This script identifies HVG by Class and also identifies TF expressed in class
suppressPackageStartupMessages({
  library(scran)
  library(scater)
})

DIR_RNA <-"~/RNA"
DIR_OUT <-"~/SCENICOut/Class"

setwd(DIR_OUT)

Cls <- commandArgs(trailingOnly = TRUE)

print(paste0("Identified HVG and TF for pyscenic run for ..", args))
## 1. --------
## obtain plotdata
load(paste0(DIR_RNA,"plotdataHBv2.RData"))
plot.data <- plot.data[plot.data$Class==args,]
## subset for atac genes
gsm <- readLines("~/extrafiles/genesGSM.txt")
pc <- read.table("~/extrafiles/gencdH38p13r37_PC.txt", header = TRUE)
pc <- as.character(pc$Gene.ID)
gsm <- gsm[gsm %in% pc]

## obtain logcounts
load(paste0(DIR_RNA,"lgcountsHB.RData"))
gsm <- gsm[gsm %in% row.names(mdt)]
mdt <- mdt[gsm,row.names(plot.data)]

## finding TFs that are expressed in at least 10% of the cells in a cluster
tfs <- readLines("~/SCENICfiles/allTFs_hg38.txt")
tfs <- tfs[tfs %in% row.names(mdt)]
mdt_tf <- mdt[tfs,]

exp <- c()
annlv2 <- unique(plot.data$Cluster)
for(x in annlv2){
  mdt_tfx <- mdt_tf[,plot.data$Cluster==x]
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

# downsample -----
tb <- data.frame(table(plot.data$Cluster))
cl <- as.character(tb$Var1)[tb$Freq>3000]
keep <- plot.data$Cluster %in% cl
table(keep)
pds1 <- plot.data[keep, ]
pds2 <- plot.data[!keep, ]
pds3 <- c()
rm(keep,tb)
##downsampling large annotations to max 3k cells
for(y in cl){
  del <- pds1[pds1$Cluster==y,]
  rn <- row.names(del)
  del2 <- del[sample(rn,3000),]
  pds3 <- rbind(pds3,del2)
  rm(del,rn,del2)
}
rm(y)
plot.data <- rbind(pds2,pds3)
mdt <- mdt[,row.names(plot.data)]

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
#hvg <- hvg[!(hvg %in% SX)] #not removing sex genes


save(hvg,heg,tfs_sel, file = paste0("~/HVGs/",Cls,"/",Cls,"_HVGnTF.RData"))
q()
