## Preparing HB normal Glial lineage + microglia data for RECORDR analysis
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")
suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(glmGamPoi)
  library(Seurat)
})
options(future.globals.maxSize = 50 * 1024^3) 

##Directories 
DIR_RNA <- "~/RNA/"
DIR_OUT <- "~/RECORDR/Normal/"

setwd(DIR_RNA)

### 1.1 load expression and metadata ------

## this is compiled raw counts expression data for the entire hindbrain
load("combinedcountsHB.RData")

## Removing ribosomal, mitochondiral and genes located on sex chromosomes
load("~/extrafles/genesnoMTRibo.RData")
load("~/extrafles/HGNCsexchr.RData")
genes <- genes[!(genes %in% SX)]

##load metadata, combined
load("plotdataHBv2.RData")
plot.data <- plot.data[plot.data$Subcluster!="ND",]
pd <- plot.data
##change dir
setwd(DIR_OUT)

# subset to glial lineage
plot.data <- plot.data[plot.data$Class %in% c("PreO","MixedGlia","Glioblasts","Immune","MOG"),]
## Remvoe ependymal cells and immune lymphocytes
cl <- unique(plot.data$Cluster)

cl_rem <-  c("Glioblasts:Gliogenic","Glioblasts:AES","Glioblasts:Progenitor_AES",
             "MixedGlia:Astrocytes_Protoplasmic","MixedGlia:Bergmann",
             "MixedGlia:Choroid_plexux","MixedGlia:Ependymal","Glioblasts:Bipotent",
             "Glioblasts:Ependoblast","Immune:Lymphocytes")


cl <- cl[!(cl %in% cl_rem)]
plot.data <- plot.data[plot.data$Cluster %in% cl,]

## remove batches with less than 200 cells
tb <- data.frame(table(plot.data$Batch))
tb <- as.character(tb$Var1)[tb$Freq>199]
plot.data <- plot.data[plot.data$Batch %in% tb,]

## I am going to create three different overlapping subset of
## clusters, to achieve more robustness in my analysis
set.seed(777)
ncells <- 3000
cls <- unique(plot.data$Cluster)
cells <- c()
for(x in cls){
  del_cells <- row.names(plot.data)[plot.data$Cluster==x]
  if(length(del_cells)>ncells){
    del_cells <- sample(del_cells,ncells)
  }
  cells <- c(cells,del_cells)
  rm(del_cells)
}
rm(x)

write(cells, file = paste0(DIR_OUT,"GlialnMicro_cells_pr_a.txt")) ##two other _b and _c
q()
cells <- readLines(paste0(DIR_OUT,"GlialnMicro_cells_pr_a.txt"))

plot.data <- plot.data[row.names(plot.data) %in% cells,]

## now remove batches with less than 100 cells
tb <- data.frame(table(plot.data$Batch))
tb <- as.character(tb$Var1)[tb$Freq>99]
plot.data <- plot.data[plot.data$Batch %in% tb,]
pd <- pd[pd$Batch %in% tb,]

## subset count matrix
mdt <- counts(com.sce)[genes,row.names(plot.data)]

## Gene filtering step, a gene needs to be expressed in at least 5% of a cells
## in a normal cell cluster

genes_sel <- NULL
genes <- row.names(mdt)
cls <- unique(plot.data$Cluster)

for(x in cls){
  del_cells <- row.names(plot.data)[plot.data$Cluster==x]
  matx_del <-mdt[,del_cells]
  rsum <- rowSums(matx_del>0)
  rsum <-(rsum/dim(matx_del)[2])*100
  genes_sel <- c(genes_sel,genes[rsum>5])
  rm(del_cells,matx_del,rsum)
}
genes_sel <- unique(genes_sel)

### susbet for protein coding genes
pc <- readLines("~/extrafiles/gencdH38p13r37_PC.txt")
genes_sel <- genes_sel[genes_sel %in% pc]
rm(mdt)
mdt <- counts(com.sce)[genes_sel,row.names(pd)]
rm(com.sce)
## using Seurat SCTtrasnform to obtain pearson residual matrix per sample
sams <- unique(pd$Batch)
R <- c()
for(x in sams){
  cell_sel <- row.names(pd)[pd$Batch==x]
  pd_del <- pd[cell_sel,]
  
  sct <- CreateSeuratObject(counts=mdt[,cell_sel],
                            meta.data = pd_del)
  sct <- SCTransform(sct,ncells = 3000, 
                     vst.flavor = "v2",
                     return.only.var.genes = FALSE)
  R_del <- sct@assays$SCT$scale.data
  gene_selx <- row.names(R_del)
  if(!all(genes_sel %in% gene_selx)){
    gene_sely <- genes_sel[!(genes_sel %in% gene_selx)]
    R_del2 <- matrix(0,length(gene_sely),dim(R_del)[2])
    row.names(R_del2) <- gene_sely
    colnames(R_del2) <- colnames(R_del)
    R_del <- rbind(R_del,R_del2)
    rm(gene_sely,R_del2)
  }
  
  R <- cbind(R,R_del[genes_sel,])
  rm(cell_sel,sct,R_del,gene_selx)
}
rm(mdt)

R <- R[genes_sel,row.names(plot.data)]

repl_python()
import numpy as np
from scipy.sparse import csr_matrix, save_npz



mat = csr_matrix(r.R)
genes_sel=r.genes_sel

save_npz(f"genemat_nor_pc_a.npz", mat) ## two other _b and _c
np.save(f"genes_pc_a.npy", genes_sel)
exit

