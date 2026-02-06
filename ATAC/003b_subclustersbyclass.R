### Subsetting ATAC by class for improved label trasnfer
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyverse)
  library(scater)
  library(scran)
})

args = commandArgs(trailingOnly=TRUE)
set.seed(456)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/ATACana/Outs/"
DIR_RNA <- "~/RNA/"
addArchRThreads(threads = 6)

setwd(DIR)
# # ## 1. load integrated ATAC to generate the Oligo subset-------
print(paste0("Creating subsetted ArchR project for class: ",args))

##not using comATACx as that is too heavy.. for Peak analysis I can do separate run
proj <- loadArchRProject("comATAC/")
load("plotdataHB_ATACv5.RData")
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$ClassX!="ND",] ##remove unassigned cells
plot.data <- plot.data.ATAC[plot.data.ATAC$ClassX==args,]
cells <- row.names(plot.data)
print("Subsetting ArchRproject by class:")

subsetArchRProject(proj,
                   cells=cells,
                   outputDirectory = paste0("Subset_",args), force = TRUE)
rm(proj, cells)
#q()


## 2. Integrate with Class RNA ------
print(paste0("Adding RNA to subsetted ArchR project for class: ",args))

proj <- loadArchRProject(paste0("Subset_",args))
load("plotdataHB_ATACv5.RData")
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$ClassX!="ND",]

cells <- getCellNames(proj)
plot.data.ATAC <- plot.data.ATAC[cells,]
proj$Clusters_subset <- plot.data.ATAC$Clusters_subset
proj$Class <- plot.data.ATAC$Class
proj$ClassX <- plot.data.ATAC$ClassX
proj$Annotationlv2_step1 <- plot.data.ATAC$Annotationlv2_step1

## loading RNA data for the class
load(paste0(DIR_RNA,"plotdataHBRNA.RData")) 


keep <- plot.data$Subset=="ND"
table(keep)
plot.data <- plot.data[!keep,]
if(args=="Neuroblasts"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("Neuroblasts_I","Neuroblasts_II","Neuroblasts_III"),]
}else if(args=="Neurons"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("Neurons_I","Neurons_II"),]
}else if(args=="AstroGlia"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("Glioblasts","MixedGlia"),]
}else if(args=="GCcom"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("GC","GCUBC"),]
}else{
  plot.data <- plot.data[plot.data$Annotationlv1==args,]
}


## loading logcounts of RNA expression 
load(paste0(DIR_RNA,"lgcountsHB.RData"))
mdt <- mdt[,row.names(plot.data)]
dec<- modelGeneVar(mdt)
dec <- data.frame(dec)
dec <- dec[!is.nan(dec$FDR),]
dec<- dec[dec$bio>0,]
dec <- dec[dec$mean>0.02,]
dec <- dec[order(dec$bio,decreasing = TRUE),]
hvg <- row.names(dec)

load(paste0(DIR_RNA,"genesnoMTRibo.RData"))
hvg <- hvg[hvg %in% genes]
load(paste0(DIR_RNA,"HGNCsexchr.RData"))
hvg <- hvg[!(hvg %in% SX)]

rm(mdt)
#preparing SE data for integration to ArchR project
gx <- readLines("extrafiles/genesGSM.txt")
gx <- gx[!duplicated(gx)]
hvg <- hvg[hvg %in% gx]
if(length(hvg)>2000){
  hvg <- hvg[1:2000]
  
}

## using raw counts
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
com.sce <- com.sce[gx,row.names(plot.data)]

##adding annotation to SCE object
com.sce$Cluster <- plot.data$Cluster ## once agan at Cluster level
rna_cells <- row.names(plot.data)

## 3. running Archr integration ---------
print("Adding gene expression:")

#add GeneIntegrationMatrix
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneExpressionMatrix",
  reducedDims = "IterativeLSI_int",
  seRNA = com.sce,
  sampleCellsATAC = 5000,##this causes division of atac datainto blocks
  sampleCellsRNA = 5000,
  dimsToUse = 1:100,
  genesUse = hvg,
  #scaleTo = 1, # 1 when using log-norm counts
  addToArrow = TRUE,
  groupRNA = "Cluster",
  nameCell = "predictedRNAlv2", #store the cell ID from the matched scRNA-seq cell
  nameGroup = "Cluster_step2", #store the group ID from the scRNA-seq cell
  nameScore = "predictedScorelv2", #store the cross-platform integration score
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
                metric = "euclidean",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = paste0("Subset_",args))
pd <- getCellColData(proj)
del <- getEmbedding(proj,"UMAP_sub")
colnames(del) <- c("iUMAP1","iUMAP2")
pd <- cbind(pd,del[row.names(pd),])

save(pd, file=paste0("Subset_",args,"/plotdata_RNAintstep2.RData"))
print("Finished!!")

q()


