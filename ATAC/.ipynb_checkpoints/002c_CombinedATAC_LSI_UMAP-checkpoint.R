## Using this script I perform LSI factorization for the integrated data
## Also call cluster and generate UMAP
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(tidyverse)
})


set.seed(456)
addArchRThreads(threads = 12)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/ATACana/Outs/"
setwd(DIR)
Sample="comATAC/"

load("plotdataHB_ATAC.RData")

################ ATAC: load the data ########################
proj <- loadArchRProject(path = Sample)

# ######################## LSI ################################
# #calculate iterative LSI
proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI_int",
                        iterations=5,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4, 0.8),
                          sampleCells = 20000, ##this causes all the cells to be used as all samples are less than 29k
                          n.start = 10
                        ),
                        varFeatures = 100000,
                        dimsToUse = 1:100,
                        totalFeatures = 500000,
                        seed = 1,
                        LSIMethod = 1,
                        scaleDims = FALSE,
                        corCutOff = 0.75,
                        excludeChr = c("chrX", "chrY", "chrMT"),
                        binarize = T,
                        force = T)

#saveArchRproj in between
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample)
# ## post LSI ------
 print("Adding clusters:")
 proj <- addClusters(input = proj,
                     name = "Clusters_int",
                     reducedDims = "IterativeLSI_int",
                     method = "Seurat",
                     force = T,
                     maxClusters = 100, #increased from 25
                     resolution=1,
                     corCutOff = 0.75,
                     scaleDims = FALSE,
                     seed = 1)


# UMAp calculation
proj <- addUMAP(ArchRProj = proj,
                name = "UMAP_int",
                reducedDims = "IterativeLSI_int",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)
#saveArchRproj
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = Sample)

print("Finished Combining Arrow Files to one proj")
plot.data.ATAC <- data.frame(getCellColData(proj))
del <- getEmbedding(proj,"UMAP_int")
colnames(del) <- c("intUMAP1","intUMAP2")
head(del)
plot.data.ATAC <- cbind(plot.data.ATAC,del)
 
plot.data.ATAC$ATAC_barcodes <- row.names(plot.data.ATAC)
save(plot.data.ATAC, file = "plotdataHB_ATACv2.RData")


suppressPackageStartupMessages({
  library(ggsci)
  library(ggrepel)
})
#load("plotdataHB_ATACv2.RData")
#define order of clusters
d <- length(unique(plot.data.ATAC$Clusters_int))
cl_lv <- paste0("C",seq(d))
plot.data.ATAC$Clusters_int<- factor(plot.data.ATAC$Clusters_int, levels = cl_lv)

##plotting cluster on UMAP----
label.d = plot.data.ATAC %>% group_by(Clusters_int) %>% 
  select(intUMAP1, intUMAP2) %>% summarize_all(median)

p <- ggplot(plot.data.ATAC, aes(x=intUMAP1, y=intUMAP2, color=Clusters_int))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Clusters_int),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


tiff(p"Figures/combined_UMAP_cluster_LSIint.tiff", unit = "in",  width=8, height=8, res = 300)
print(p)
dev.off()
rm(p)


##plotting Sample on UMAP----
d <- length(unique(plot.data.ATAC$Sample))
label.d = plot.data.ATAC %>% group_by(Sample) %>% 
  select(intUMAP1, intUMAP2) %>% summarize_all(median)
#for plot
p <- ggplot(plot.data.ATAC, aes(x=intUMAP1, y=intUMAP2, color=Sample))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Sample),size = 2.5, data = label.d, show.legend = FALSE, max.overlaps = 25)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


tiff("Figures/combined_UMAP_cluster_Sample.tiff", unit = "in",  width=8, height=8, res = 300)
print(p)
dev.off()
rm(p)



