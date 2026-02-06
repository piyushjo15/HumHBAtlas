### Subsetting ATAC clusters by groups and adding gene expression data
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/.conda/envs/pjexv4/bin/python3.7")
use_python("~/.conda/envs/pjexv4/bin/python3.7")

suppressPackageStartupMessages({
  library(ArchR)
  library(hexbin)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(dplyr)
  library(tidyr)
})

args = commandArgs(trailingOnly=TRUE)

set.seed(456)
addArchRThreads(threads = 12)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "/home/p541i/MBsnANA/HBana/ATACana/Outs/"
DIRS <- "/home/p541i/MBsnANA/HBana/ATACana/HBATAC_Scripts/"
DIR_RNA <- "/home/p541i/DATA/Diemana/scepro/"

setwd(DIR)
Sample="comATAC/"


############### ATAC: load the data ########################
proj <- loadArchRProject(path = Sample)

###sub-setting by  clusters associated with celltype
load(paste0(DIR,"extrafiles/plotdataHB_ATACv4.RData"))
cl_lv <- readLines(paste0(DIRS,args,".txt"))
keep <- plot.data.ATAC$Clusters_subset %in% cl_lv
table(keep)
plot.data.ATAC <- plot.data.ATAC[keep,]
atac_cells <- row.names(plot.data.ATAC)

print("Subsetting ArchRproject by supercluster:")

subsetArchRProject(proj,
                   cells=atac_cells,
                   outputDirectory = paste0("Subset_",args))
rm(proj)

####### Adding gene expression data ############
##load subsetted archr project
proj <- loadArchRProject(paste0("Subset_",args))
proj <- proj[atac_cells,] ##making sure order is correct
proj$Clusters_subset <- plot.data.ATAC$Clusters_subset

##loading counts and  logcounts gene expression
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
load(paste0(DIR_RNA,"plotdata",args,"com.RData")) #add com for GC
##remove ND
plot.data <- plot.data[plot.data$Annotationlv3!="ND",]
#cts <- round(counts(com.sce[,row.names(plot.data)]))
#rm(com.sce)
# load(paste0(DIR_RNA,"lgcountsHB.RData"))
# mdt <- mdt[,row.names(plot.data)]

#preparing SE data for integration to ArchR project
##I got error in Seurat Findanchors, so I now I am subsetting to only genes present in GSM matrix
gx <- readLines("extrafiles/genesGSM.txt")
gx <- gx[!duplicated(gx)]
com.sce <- com.sce[gx,row.names(plot.data)]
com.sce$Annotationlv3 <- plot.data$Annotationlv3
rna_cells <- row.names(plot.data)
#rm(cts,mdt)
load(paste0(DIR_RNA,"GCcomhvg.RData")) #Oligocom_psb.RData for Oligo
rm(del)
hvg <- hvg[hvg %in% gx]
print("Adding gene expression:")

#add GeneIntegrationMatrix
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneExpressionMatrix",
  reducedDims = "IterativeLSI_int",
  seRNA = com.sce,
  sampleCellsATAC = 15000,##this causes division of atac datainto blocks
  sampleCellsRNA = 20000,
  dimsToUse = 1:50,
  genesUse = hvg,
  addToArrow = TRUE,
  groupRNA = "Annotationlv3",
  nameCell = "predicted_RNA", #store the cell ID from the matched scRNA-seq cell
  nameGroup = "predictedGroup", #store the group ID from the scRNA-seq cell
  nameScore = "predictedScore", #store the cross-platform integration score
  force = TRUE ,
  useImputation = FALSE
)


##adding imputed weights
print("Adding imputed weights:")

proj <- addImputeWeights(proj, reducedDims = "IterativeLSI_int")

#
print("Adding UMAP:")

proj <- addUMAP(ArchRProj = proj,
                name = "UMAP_sub",
                reducedDims = "IterativeLSI_int",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = paste0("Subset_",args))
# post saving -----------
proj <- loadArchRProject(paste0("Subset_",args))

# ##extract LSI and store as .npz
# LSI <- getReducedDims(proj,"IterativeLSI_int")
# cn <- row.names(LSI)
# cn <- gsub("-1","",cn)
# cn <- gsub("#","_",cn)
# row.names(LSI) <- cn
# repl_python()
# import numpy as np
# import scipy.sparse as sp
# import os
# Sample=r.args
# projDir = os.path.join("/home/p541i/MBsnANA/HBana/ATACana/Outs/Subset_"+Sample+"/")
# LSI=sp.csr_matrix(r.LSI)
# sp.save_npz(os.path.join(projDir,"LSImat.npz"), LSI)
# exit

##extracting GSM for subset-
#mdt <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
#save(mdt, file = paste0("Subset_",args,"/GSMmat.RData"))
#q()
##fixing labels 
plot.data.ATAC <- data.frame(getCellColData(proj))
cM <- confusionMatrix(plot.data.ATAC$Clusters_subset, plot.data.ATAC$predictedGroup)
labelOld <- rownames(cM)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
df <- data.frame(Old=labelOld, New=labelNew)
row.names(df) <- df$Old
plot.data.ATAC$predictedGroup_fixed <- "ND"
cls <- unique(plot.data.ATAC$Clusters_subset)
for(x in cls){
  keep <- plot.data.ATAC$Clusters_subset==x
  plot.data.ATAC[keep,"predictedGroup_fixed"] <- df[x,"New"]
  rm(keep)
}
rm(x)
table(plot.data.ATAC$predictedGroup_fixed)
del <- getEmbedding(proj,"UMAP_sub")
colnames(del) <- c("UMAP1","UMAP2")
head(del)
plot.data.ATAC <- cbind(plot.data.ATAC,del)
save(plot.data.ATAC, file = paste0("extrafiles/plotdataATAC_",args,".RData"))

q()
load("extrafiles/plotdataATAC_GC.RData")
table(plot.data.ATAC$predictedGroup)
