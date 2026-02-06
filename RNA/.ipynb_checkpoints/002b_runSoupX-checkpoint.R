## This script used cell identified from DIEM and performs SoupX and decontX based
## debris signal correciton

suppressPackageStartupMessages({
  library(diem)
  library(ggplot2)
  library(gridExtra)
  library(Matrix)
  library(scater)
  library(scran)
})
#argument
args = commandArgs(trailingOnly=TRUE)
#Load the premrna counts
DIR="rawmatrix/"
counts.in <- read_10x(paste0(DIR,args[1],"/GeneFull"))
dim(counts.in)
# load the exon part

counts.ex <- read_10x(paste0(DIR,args[1],"/Gene"))
dim(counts.ex)
#table(row.names(counts.in)==row.names(counts.ex))
#check colnames and rownames
if(all(colnames(counts.in) ==colnames(counts.ex)) & all(row.names(counts.in) == row.names(counts.ex))){
  
  #head(counts.ex)[,1:10]
  #head(counts.ex)[,1:10]
  genes <- read.delim("~/HB/ann/gencdH38p13r37CR_genesN.txt", header = FALSE)
  genes <- unique(genes$V1)
  
  #finding genes that were counted less with intronic marix
  mean.in <- rowMeans(counts.in)
  mean.ex <- rowMeans(counts.ex)
  d <- (mean.ex > mean.in)
  print(table(d))
  #replacing above genes' data in intronic matrix from exonic matrix
  counts.in2 <- counts.in[!d,]
  dim(counts.in2)
  counts.ex2 <- counts.ex[d,]
  dim(counts.ex2)
  
  counts.in2 <- rbind(counts.in2,counts.ex2)
  dim(counts.in2)
  counts.in2 <- counts.in2[genes,]
  rm(genes)
}


rm(counts.ex2, counts.ex, counts.in)
load(paste0("diemv4",args[1],".RData")) ##output of 002a_rundiem.R
cells <- row.names(mdt.pd[mdt.pd$Call=="Clean",])
sce <- SingleCellExperiment(list(counts=mdt@counts[,cells]))
#table(row.names(drop_data)==colnames(sce))
del <- mdt.pd[colnames(sce),]
colData(sce) <- DataFrame(del)

mdfr3 <- mdfr2[colnames(sce),]
sce$in_ex_frac <- mdfr3$frac

row.names(mdt.pd) <- paste0(args[1],"_",row.names(mdt.pd))
colnames(sce) <- paste0(args[1],"_",colnames(sce))

#head(colData(sce))
#initial filtering based on pct.mt and feature counts,
# already filtered based on min_gene =500 in diem run
x <- mean(sce$pct.mt)
z <- x+ 1*mad(sce$pct.mt)
if (z <2) {z=2}
##based on boxplot distribution
co <- quantile(sce$nUMIs, 0.975)
keep <- (sce$pct.mt < z) & (sce$nUMIs < co)
sce <- sce[,keep]
keep <- sce$in_ex_frac >0.25
sce <- sce[,keep]

rm(x,keep,z, mdt, mdfr3)

#removing soup----
suppressPackageStartupMessages({
  library(SoupX)
  library(Seurat)
  library(dplyr)
})

#to select important debris genes
#find custers and UMAP
sce.S <- CreateSeuratObject(counts = counts(sce))
sce.S <- SCTransform(sce.S, verbose = FALSE)
sce.S <- RunPCA(sce.S)
sce.S <- RunUMAP(sce.S, dims = 1:30)
sce.S <- FindNeighbors(sce.S, dims = 1:30)
sce.S <- FindClusters(sce.S, resolution = 1)

filt.matrix <- counts(sce)
soup.channel  <- SoupChannel(counts.in2, filt.matrix)
rm(counts.in2)
meta    <- sce.S@meta.data
umap    <- sce.S@reductions$umap@cell.embeddings
soup.channel  <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
soup.channel  <- setDR(soup.channel, umap)
#head(meta)
soup.channel  <- autoEstCont(soup.channel, tfidfMin = 0.8, soupQuantile = 0.85, forceAccept = TRUE)

soup.data <- soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ]

soup.data$symbol <- row.names(soup.data)

tiff(paste0(args[1],"_SoupX_geneplot.tiff"),width = 6, height = 4, units ="in", res=150)

ggplot(filter(soup.data, rank(-est) < 100) , aes(x=rank(-est), y=est, label=symbol)) +
  geom_point(size=0.1, alpha=0.3) +
  geom_text(check_overlap = T) +
  theme_bw()

dev.off()


soupx_corrected  <- adjustCounts(soup.channel, roundToInt = T)
assay(sce,"unfiltcounts") <- counts(sce)
counts(sce) <- soupx_corrected

library(celda)
#deconTX
sce <- decontX(sce, z = sce.S$seurat_clusters)
assay(sce,"rawcounts") <- counts(sce)
counts(sce) <-  decontXcounts(sce)
decontXcounts(sce) <- NULL

save(sce,mdfr2,mdt.pd, file = paste0(args[1],"_postSoupXDX.RData"))

q()
