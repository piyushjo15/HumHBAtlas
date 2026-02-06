## Label transfer round 1 from integrated RNA to integrated ATAC
## 
suppressPackageStartupMessages({
  library(ArchR)
  library(scran)
  library(scater)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(tidyverse)
})

set.seed(456)
addArchRThreads(threads = 6)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/ATACana/Outs/"
DIR_RNA <- "~/RNA/"

setwd(DIR)
Sample="comATACx/" ### i copied 'comATAC' into 'comATACx' to be safe and used this for further analysis

proj <- loadArchRProject(path = Sample)
## 1. ---------------
##loading counts and  logcounts gene expression
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
load(paste0(DIR_RNA,"plotdataHBRNA.RData")) 
keep <- plot.data$Subcluster=="ND" ##removed cells not annotated
table(keep)
plot.data <- plot.data[!keep,]
#preparing SE data for integration to ArchR project
##I got error in Seurat Findanchors, so I now I am subsetting to only genes present in GSM matrix
gx <- readLines("extrafiles/genesGSM.txt") ##These are geenes present in ATAC gene-score matrix
gx <- gx[!duplicated(gx)]
com.sce <- com.sce[gx,row.names(plot.data)]
com.sce$Cluster <- plot.data$Cluster ##going to trasnfer labels at cluster level
rna_cells <- row.names(plot.data)
load(paste0(DIR_RNA,"topHVGcomhvgnoRiboMTSX.RData"))  ##HVG from normal cells

hvg <- topHVG[topHVG%in% gx] ##subset to ATAC genes
print("Adding gene expression:")


#add GeneIntegrationMatrix
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneExpressionMatrix",
  reducedDims = "IterativeLSI_int",
  seRNA = com.sce,
  sampleCellsATAC = 20000,##this causes division of atac datainto blocks
  sampleCellsRNA = 30000,
  dimsToUse = 1:100,
  genesUse = hvg[1:3000],
  addToArrow = FALSE,
  groupRNA = "Cluster",
  nameCell = "predRNA", #store the cell ID from the matched scRNA-seq cell
  nameGroup = "predCluster", #store the group ID from the scRNA-seq cell
  nameScore = "predScore", #store the cross-platform integration score
  force = TRUE ,
  useImputation = FALSE
)

print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory =Sample)
#q()
## 2. --------
##fixing labels 
load("plotdataHB_ATACv3.RData") ##this is plot.data.ATAC
plot.data <- data.frame(getCellColData(proj))
del <- getEmbedding(proj,"UMAP_int")
colnames(del) <- c("UMAP1","UMAP2")
plot.data <- cbind(plot.data,del[row.names(plot.data),])
table(row.names(plot.data) %in% row.names(plot.data.ATAC))
plot.data.ATAC <- plot.data.ATAC[row.names(plot.data),]
plot.data.ATAC$Cluster_step1 <- plot.data$predCluster
plot.data.ATAC$UMAP1 <- plot.data$UMAP1
plot.data.ATAC$UMAP2 <- plot.data$UMAP2
plot.data.ATAC$predictedGroup <- 
  plot.data.ATAC$predictedScore <- NULL

## 'Cluster_subset' are ATAC clustered identified at lv2 clustering
## here I am trying to identify which ATAC cluster matches to which RNA cluster
cM <- confusionMatrix(plot.data.ATAC$Clusters_subset, plot.data.ATAC$Cluster_step1) 
labelOld <- rownames(cM)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
df <- data.frame(Old=labelOld, New=labelNew)
row.names(df) <- df$Old
write.table(df,"extrafiles/Clusters_subset_lv2.txt", sep = "\t",row.names = FALSE,
            quote = FALSE)
## get class from rna
sub2cls <- read.delim(paste0(DIR_RNA,"Classcluster.txt")) ## Assgining RNA Class to ATAC clusters
row.names(sub2cls) <- sub2cls$Cluster
head(sub2cls)
plot.data.ATAC$Class_pred <- "ND"
clx <- unique(plot.data.ATAC$Cluster_step1)
for(x in clx){
  plot.data.ATAC[plot.data.ATAC$Cluster_step1==x,"Class_pred"] <- sub2cls[x,"Class"]
}
## Also doing a sanity check to identify dominant class per ATAC cluster
## This doesn't automatically become Class assigned to ATAC cells
## but is just needed for comparison
cM <- confusionMatrix(plot.data.ATAC$Clusters_subset, plot.data.ATAC$Class_pred)
labelOld <- rownames(cM)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
df <- data.frame(Old=labelOld, New=labelNew)
row.names(df) <- df$Old
write.table(df,"extrafiles/Clusters_subset_Class.txt", sep = "\t",row.names = FALSE,
            quote = FALSE)
save(plot.data.ATAC, file = "plotdataHB_ATACv4.RData")
q()
### Here I verfied the Class that should be assigned to each of the ATAC clauster 
## fix class, after mannual inspection
load("plotdataHB_ATACv4.RData") ## this is plot.data.ATAC
ann <- read.delim("Cluster_sub_ANN.txt")
row.names(ann) <- ann$Clusters_subset
head(ann)
plot.data.ATAC$Class <- "ND"
clx <- unique(plot.data.ATAC$Clusters_subset)
for(x in clx){
  plot.data.ATAC[plot.data.ATAC$Clusters_subset==x,"Class"] <- ann[x,"Class"]
}
df <- data.frame(table(plot.data.ATAC$Clusters_subset,plot.data.ATAC$Class))
df <- df[df$Freq>0,]
## for integration at step 2, maybe mixing GB+MixedGlia as one, all neuroblasts
## as NB and all Neurons as NC is better?
plot.data.ATAC$ClassX <- plot.data.ATAC$Class

keep <- plot.data.ATAC$Class %in% c("Glioblasts","MixedGlia")
plot.data.ATAC[keep,"ClassX"] <- "AstroGlia"

keep <- plot.data.ATAC$Class %in% c("Neuroblasts_I","Neuroblasts_II","Neuroblasts_III")
plot.data.ATAC[keep,"ClassX"] <- "Neuroblasts"

keep <- plot.data.ATAC$Class %in% c("Neurons_I","Neurons_II")
plot.data.ATAC[keep,"ClassX"] <- "Neurons"

keep <- plot.data.ATAC$Class %in% c("GCUBC","GC")
plot.data.ATAC[keep,"ClassX"] <- "GCcom"

save(plot.data.ATAC, file = "plotdataHB_ATACv5.RData")

