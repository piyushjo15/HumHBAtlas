### Subsetting ATAC by class for improved analysis
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
  library(tidyr)
  library(scater)
  library(scran)
})

args = commandArgs(trailingOnly=TRUE)
set.seed(456)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/MBsnANA/HBana/ATACana/Outs/"
DIRS <- "~/MBsnANA/HBana/ATACana/HBATAC_Scripts/"
DIR_RNA <- "~/MBsnANA/Diemana/scepro/"
addArchRThreads(threads = 6)

setwd(DIR)

# # # ## 1. load integrated ATAC to generate the Oligo subset-------
# print(paste0("Creating subsetted ArchR project for class: ",args))
# 
# ##not using comATACx as that is too heavy.. for Peak analysis I can do separate run
# proj <- loadArchRProject("comATAC/")
# load("extrafiles/plotdataHB_ATACv6_260325.RData")
# plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$ClassX!="ND",]
# plot.data <- plot.data.ATAC[plot.data.ATAC$ClassX==args,]
# 
# # #remove samples less than 10 cells,
# # # doing this for  MES_VAS
# # tb <- data.frame(table(plot.data$Sample))
# # tb <- as.character(tb$Var1)[tb$Freq>9]
# # plot.data <- plot.data[plot.data$Sample %in% tb,]
# 
# cells <- row.names(plot.data)
# print("Subsetting ArchRproject by class:")
# 
# subsetArchRProject(proj,
#                    cells=cells,
#                    outputDirectory = paste0("Subset_",args), force = TRUE)
# rm(proj, cells)
# q()
#
## 2. Integrate with Class RNA ------
print(paste0("Adding RNA to subsetted ArchR project for class: ",args))

proj <- loadArchRProject(paste0("Subset_",args))
load("extrafiles/plotdataHB_ATACv6_260325.RData")
plot.data.ATAC <- plot.data.ATAC[plot.data.ATAC$ClassX!="ND",]

cells <- getCellNames(proj)
plot.data.ATAC <- plot.data.ATAC[cells,]
proj$Clusters_subset <- plot.data.ATAC$Clusters_subset
proj$Class <- plot.data.ATAC$Class
proj$ClassX <- plot.data.ATAC$ClassX
proj$Annotationlv2_step1 <- plot.data.ATAC$Annotationlv2_step1

## loading RNA data for the class
load(paste0(DIR_RNA,"plotdataHBTempANN250305.RData")) 


keep <- plot.data$Annotationlv3=="ND"
table(keep)
plot.data <- plot.data[!keep,]
if(args=="Neuroblasts"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("Neuroblasts_A","Neuroblasts_B","Neuroblasts_C"),]
}else if(args=="Neurons"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("Neurons_A","Neurons_B"),]
}else if(args=="AstroGlia"){
  plot.data <- plot.data[plot.data$Annotationlv1 %in% c("Glioblast","MixedGlia"),]
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

# #using lognorm counts
# ##both raw counts and lognorm with scale=1 look very similar
# mdt <- mdt[gx,row.names(plot.data)]
# com.sce <- SingleCellExperiment(list(logcounts=mdt))

##adding annotation to SCE object
com.sce$Annotationlv2 <- plot.data$Annotationlv2
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
  groupRNA = "Annotationlv2",
  nameCell = "predictedRNAlv2", #store the cell ID from the matched scRNA-seq cell
  nameGroup = "Annotationlv2_step2", #store the group ID from the scRNA-seq cell
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
## 4. plotting---------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggsci)
})
args <- "Neuroblasts"
load(paste0("Subset_",args,"/plotdata_RNAintstep2.RData"))
pd <- data.frame(pd)

d <- length(unique(pd$Annotationlv2_step1))
label.d = pd %>% group_by(Annotationlv2_step1) %>%
  select(iUMAP1, iUMAP2) %>% summarize_all(median)

# #for plot
p1 <- ggplot(pd, aes(x=iUMAP1, y=iUMAP2, color=Annotationlv2_step1))+
  geom_point(size=1)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Annotationlv2_step1),size = 2.5, 
                   data = label.d, show.legend = FALSE, max.overlaps = 50)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")
d <- length(unique(pd$Annotationlv2_step2))
label.d = pd %>% group_by(Annotationlv2_step2) %>%
  select(iUMAP1, iUMAP2) %>% summarize_all(median)

# #for plot
p2 <- ggplot(pd, aes(x=iUMAP1, y=iUMAP2, color=Annotationlv2_step2))+
  geom_point(size=1)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Annotationlv2_step2),size = 2.5, 
                   data = label.d, show.legend = FALSE, max.overlaps = 50)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")
gridExtra::grid.arrange(p1,p2)
sort(table(pd$Annotationlv2_step2),decreasing = TRUE)

##when step 2 annotation transfer is done at lv3
# ##convert Annlv3 to Annlv2
# sub2clst <- read.delim("~/MBsnANA/Diemana/scepro/Classcluster.txt")
# row.names(sub2clst) <- sub2clst$Annlv3
# pda$predlv2 <- "ND"
# clx <- unique(pda$predictedGrouplv3)
# for(x in clx){
#   pda[pda$predictedGrouplv3==x,"predlv2"] <- sub2clst[x,"AnnLv2new"]
# }
# pda <- pda[,c("Sample","predictedScorelv3","predlv2","UMAP1","UMAP2","CellID")]
# head(pda)
# save(pd, file=paste0("Subset_",args,"/pred_",args,"_seurat_pd.RData"))


