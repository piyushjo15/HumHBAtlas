##addinging gene integration to entire integrated ATAC
suppressPackageStartupMessages({
  library(ArchR)
  library(scran)
  library(scater)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(dplyr)
  library(tidyr)
})

set.seed(456)
addArchRThreads(threads = 6)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/MBsnANA/HBana/ATACana/Outs/"
DIRS <- "~/MBsnANA/HBana/ATACana/HBATAC_Scripts/"
DIR_RNA <- "~/MBsnANA/Diemana/scepro/"

setwd(DIR)
Sample="comATACx/"

proj <- loadArchRProject(path = Sample)
## 1. ---------------
##loading counts and  logcounts gene expression
load(paste0(DIR_RNA,"combinedcountsHB.RData"))
load(paste0(DIR_RNA,"plotdataHBTempANN250305.RData")) 
keep <- plot.data$Annotationlv3=="ND"
table(keep)
plot.data <- plot.data[!keep,]
#preparing SE data for integration to ArchR project
##I got error in Seurat Findanchors, so I now I am subsetting to only genes present in GSM matrix
gx <- readLines("extrafiles/genesGSM.txt")
gx <- gx[!duplicated(gx)]
com.sce <- com.sce[gx,row.names(plot.data)]
com.sce$Annotationlv2 <- plot.data$Annotationlv2
rna_cells <- row.names(plot.data)
load(paste0(DIR_RNA,"topHVGcomhvgnoRiboMTSX.RData")) 

hvg <- topHVG[topHVG%in% gx]
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
  groupRNA = "Annotationlv2",
  nameCell = "predRNA", #store the cell ID from the matched scRNA-seq cell
  nameGroup = "predAnnotationlv2", #store the group ID from the scRNA-seq cell
  nameScore = "predScore", #store the cross-platform integration score
  force = TRUE ,
  useImputation = FALSE
)

print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory =Sample)
q()
## 2. --------
##fixing labels 
load("extrafiles/plotdataHB_ATACv5.RData")
proj <- loadArchRProject("comATACx/")
plot.data <- data.frame(getCellColData(proj))
del <- getEmbedding(proj,"UMAP_int")
colnames(del) <- c("UMAP1","UMAP2")
plot.data <- cbind(plot.data,del[row.names(plot.data),])
table(row.names(plot.data) %in% row.names(plot.data.ATAC))
plot.data.ATAC <- plot.data.ATAC[row.names(plot.data),]
plot.data.ATAC$Annotationlv2_step1 <- plot.data$predAnnotationlv2
plot.data.ATAC$UMAP1 <- plot.data$UMAP1
plot.data.ATAC$UMAP2 <- plot.data$UMAP2
plot.data.ATAC$Annotationx <- plot.data.ATAC$predictedGroup <- 
  plot.data.ATAC$predictedScore <- plot.data.ATAC$predictedGroup_fixed <- NULL

cM <- confusionMatrix(plot.data.ATAC$Clusters_subset, plot.data.ATAC$Annotationlv2_step1)
labelOld <- rownames(cM)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
df <- data.frame(Old=labelOld, New=labelNew)
row.names(df) <- df$Old
write.table(df,"extrafiles/Clusters_subset_Annlv2.txt", sep = "\t",row.names = FALSE,
            quote = FALSE)
## get class from rna
sub2cls <- read.delim(paste0(DIR_RNA,"Classcluster.txt"))
row.names(sub2cls) <- sub2cls$Cluster
head(sub2cls)
plot.data.ATAC$Class_pred <- "ND"
clx <- unique(plot.data.ATAC$Annotationlv2_step1)
for(x in clx){
  plot.data.ATAC[plot.data.ATAC$Annotationlv2_step1==x,"Class_pred"] <- sub2cls[x,"Class"]
}
cM <- confusionMatrix(plot.data.ATAC$Clusters_subset, plot.data.ATAC$Class_pred)
labelOld <- rownames(cM)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
df <- data.frame(Old=labelOld, New=labelNew)
row.names(df) <- df$Old
write.table(df,"extrafiles/Clusters_subset_Class.txt", sep = "\t",row.names = FALSE,
            quote = FALSE)
save(plot.data.ATAC, file = "extrafiles/plotdataHB_ATACv6_260325.RData")
load("extrafiles/plotdataHB_ATACv6_260325.RData")

## fix class, after mannual inspection
ann <- read.delim("extrafiles/Cluster_sub_ANN.txt")
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

keep <- plot.data.ATAC$Class %in% c("Glioblast","MixedGlia")
plot.data.ATAC[keep,"ClassX"] <- "AstroGlia"

keep <- plot.data.ATAC$Class %in% c("Neuroblasts_A","Neuroblasts_B","Neuroblasts_C")
plot.data.ATAC[keep,"ClassX"] <- "Neuroblasts"

keep <- plot.data.ATAC$Class %in% c("Neurons_A","Neurons_B","Neurons_C")
plot.data.ATAC[keep,"ClassX"] <- "Neurons"

save(plot.data.ATAC, file = "extrafiles/plotdataHB_ATACv6_260325.RData")

load("extrafiles/plotdataHB_ATACv6_260325.RData")

# load("extrafiles/plotdataHB_ATACv6_240125.RData")
# df <-data.frame(table(plot.data.ATAC$Annotationlv2))
# df2 <-data.frame(table(plot.data.ATAC$Annotationlv2_fixed))
# row.names(df) <- as.character(df$Var1)
# row.names(df2)  <- as.character(df2$Var1)
# df$Freq2 <- 0
# df[row.names(df2),"Freq2"] <- df2$Freq
# head(df2)

d <- length(unique(plot.data.ATAC$Class))
#plot.data$NND_cl <- factor(plot.data$NND_cl, levels = paste0("NNDJS_",c(0:(d-1))))
label.d = plot.data.ATAC %>% group_by(predictedGroup) %>% 
  select(UMAP1, UMAP2) %>% summarize_all(median)
#for plot
p <- ggplot(plot.data.ATAC, aes(x=UMAP1, y=UMAP2, color=predictedGroup))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_igv()(51))(d))+
  theme_classic()+
  geom_label_repel(aes(label = predictedGroup),size = 2.5, data = label.d, show.legend = FALSE, max.overlaps = 50)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")
tiff("predRNA.tiff", unit = "in",  width=20, height=16, res = 150)
print(p)
dev.off()

## new calss on 110225

