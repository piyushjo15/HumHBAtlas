
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})


sample = commandArgs(trailingOnly=TRUE)
set.seed(456)
addArchRThreads(threads =8) 
quant_cut <- function(x){
  min_q <-quantile(x,0.01)
  max_q <-quantile(x,0.99)
  x[x<min_q] <- min_q
  x[x>max_q] <- max_q
  return(x)
}

##Dirs
DIR_ATACana <- "~/ATACana/"
DIR_FQ <- "~/snATACana/"
DIR_OUT <- "~/ATACana/Outs/"

#Next, I load the reference genome
addArchRGenome("hg38")
setwd(paste0(DIR_OUT))

# Removing doublets------
# For this we will be using simulated doublets and identifying their neighbours in the LSI space.
# Throughout this project, we are using the TF-(logIDF) LSI method as introduced by Cusanovich et al.
# We don't scale LSI components (i.e. each is proportional to its variance) but remove components with more than 0.75 correlation with sequencing depth.
# Here, we will be using 50 dimensions, 10 iterations and label the 10 NNs to each doublet.
##already created Arrow file before
print(paste0("processing sample ",sample," for ArchR pre-processing step!"))
ArrowFiles <- paste0(DIR_OUT,"Arrows/",sample,".arrow")

print("calculating doublet scores .....")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  dimsToUse = 1:50,
  nTrials = 10,
  scaleDims = F,
  threads = 1,
  LSIParams = list(seed = 1,
  varFeatures = 50000,
  excludeChr = c("chrX", "chrY", "chrMT")),
  outDir = "002_doublets"
)
#Finally, an ArchR project can be generated.
print("creating ArchR project post doublet scores .....")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = paste0("projs/",sample), 	# generate an output directory before with this name
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
# ##loading meta data, will add later
# met_sams <- read.delim("met_sams.txt",row.names = 1)
# 
# proj$Stage <- met_sams[sample,"Stage"]
# proj$Age <- met_sams[sample,"Age"]
# proj$Sex <- met_sams[sample,"Sex"]

##number of cells post intial QC
print("Number of cells before doublet filtering:")
nCells(proj)
plot.data <- data.frame(getCellColData(proj))
plot.data$barcode <- row.names(plot.data)
# Filtering putative doublets (pass 1).
# Here we want to be relatively lenient, so we select a ratio of 1 (i.e. filtering top 5% of barcodes for a sample of 5,000 cells).
# We are using Doublet enrichment as a cutoff.

proj <- filterDoublets(ArchRProj = proj,filterRatio = 1)

print("Number of cells post doublet filtering:")
nCells(proj)
##saving doublet cells
doublets <- plot.data$barcode[!(plot.data$barcode %in% getCellNames(proj))]
write(doublets, paste0("002_doublets/",sample,"/",sample,"_doublets.txt"))

z <- quantile(plot.data$nFrags,0.99) ##based on quantile
keep1 <- plot.data$nFrags > z
#remove doublet enrichment
keep2 <- plot.data$DoubletEnrichment > 20 
put_doublets <- plot.data$barcode[keep1 | keep2]
length(put_doublets)


#total putative doublets that will be removed
cells_rem <- unique(c(put_doublets,doublets))
write(cells_rem,  paste0("002_doublets/",sample,"/",sample,"_rem_cells.txt"))
keep <- plot.data$barcode %in% cells_rem
plot.data <- plot.data[!keep,]
##subsetting object to filtered cells
proj <- proj[plot.data$barcode,]


#Now we can calculate our first iterative LSI. We are now extending our dimensions to 100 and the number of variable features to 100K.

proj <- addIterativeLSI(ArchRProj = proj,
                        useMatrix = "TileMatrix",
                        name = "IterativeLSI_qc",
                        iterations=5,
                        clusterParams = list(
                          resolution = c(0.1, 0.2, 0.4, 0.8),
                          sampleCells = 20000, ##this causes all the cells to be used as all samples are less than 29k
                          n.start = 10),
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


### post LSI ------
#proj <- loadArchRProject(sample)
##iidentifying clusters to perform second QC pass with high resolution=2
print("Adding clusters:")

proj <- addClusters(input = proj,
    name = "Clusters_qc",
    reducedDims = "IterativeLSI_qc",
    method = "Seurat",
    force = T,
    resolution=2,
    corCutOff = 0.75,
    scaleDims = FALSE,
    seed = 1)

# UMAp calculation

proj <- addUMAP(ArchRProj = proj,
                name = "UMAP_qc",
                reducedDims = "IterativeLSI_qc",
                minDist = 0.2,
                metric = "cosine",
                nNeighbors = 25,
                force = T,
                seed = 1,
                scaleDims = F,
                corCutOff = 0.75)


print("saving archr project")
saveArchRProject(proj)

##extract LSI
LSI <- getReducedDims(proj,"IterativeLSI_qc")

save(LSI, file = paste0("projs/",sample,"/",sample,"_LSI.RData"))

#q()

##extract gene score matrix
mdt <- getMatrixFromProject(proj, "GeneScoreMatrix")
genx <- data.frame(rowData(mdt))
mdt <- assay(mdt)
row.names(mdt) <- genx$name

save(mdt, file = paste0("projs/",sample,"/",sample,"_GSM.RData"))


del <- getEmbedding(proj,"UMAP_qc")
colnames(del) <- c("UMAP1","UMAP2")
head(del)
plot.data <- data.frame(getCellColData(proj))
plot.data <- cbind(plot.data, del)
##plotting cluster on UMAP-----
d <- length(unique(plot.data$Clusters_qc))

label.d = plot.data %>% group_by(Clusters_qc) %>%
  select(UMAP1, UMAP2) %>% summarize_all(median)
#for plot
p <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=Clusters_qc))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Clusters_qc),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")


tiff(paste0(DIR_OUT,"Figures/",sample,"/",sample,"_UMAP_cluster.tiff"), unit = "in",  width=4, height=4, res = 150)
print(p)
dev.off()
rm(p)

## plot QC UMAP-----
#quantile trasnformation
plot.data$XX <- quant_cut(plot.data$TSSEnrichment)
p <- ggplot(plot.data%>%arrange(XX), aes(x=UMAP1, y=UMAP2, color=XX))+
  geom_point(size=0.5)+
  scale_color_viridis_c(direction = -1,option = "A", name="TSSEnrichment")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank())
tiff(paste0(DIR_OUT,"Figures/",sample,"/",sample,"_TSSEnrichment.tiff"), unit = "in",  width=5, height=4, res = 150)
print(p)
dev.off()


## Doublet score
plot.data$XX <- quant_cut(plot.data$DoubletScore)
p <- ggplot(plot.data%>%arrange(XX), aes(x=UMAP1, y=UMAP2, color=XX))+
  geom_point(size=0.5)+
  scale_color_viridis_c(direction = -1,option = "A", name="DoubletScore")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank())
tiff(paste0(DIR_OUT,"Figures/",sample,"/",sample,"_DoubletScores.tiff"), unit = "in",  width=5, height=4, res = 150)
print(p)
dev.off()

## Doublet enrichment
plot.data$XX <- quant_cut(plot.data$DoubletEnrichment)
p <- ggplot(plot.data%>%arrange(XX), aes(x=UMAP1, y=UMAP2, color=XX))+
  geom_point(size=0.5)+
  scale_color_viridis_c(direction = -1,option = "A", name="DoubletEnrichment")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank())
tiff(paste0(DIR_OUT,"Figures/",sample,"/",sample,"_DoubletEnrichment.tiff"), unit = "in",  width=5, height=4, res = 150)
print(p)
dev.off()

## nFragments
plot.data$XX <- quant_cut(plot.data$nFrags)
p <- ggplot(plot.data%>%arrange(XX), aes(x=UMAP1, y=UMAP2, color=XX))+
  geom_point(size=0.5)+
  scale_color_viridis_c(direction = -1,option = "A", name="nFrags")+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank())
tiff(paste0(DIR_OUT,"Figures/",sample,"/",sample,"_nFrags.tiff"), unit = "in",  width=5, height=4, res = 150)
print(p)
dev.off()


d <- length(unique(plot.data$Clusters_qc))
p <- ggplot(plot.data, aes(x=Clusters_qc, y=log10(nFrags), fill=Clusters_qc)) +
  geom_violin() +
  geom_boxplot(notch = T, width=0.2, alpha=0) +
  scale_fill_manual(values =colorRampPalette(pal_simpsons()(16))(d)) +
  theme_classic()+
  ylab("log10 number of fragments") +
  xlab("Clusters") +
  theme(axis.text.x = element_text(angle=45, hjust = 0.9))

tiff(paste0(DIR_OUT,"Figures/",sample,"/",sample,"_SpuriousClusters_nFrag_distribution.tiff"), unit = "in",  width=8, height=6, res = 150)
print(p)
dev.off()
save(plot.data, file =paste0(DIR_OUT,"projs/",sample,"/",sample,"_barcode_stats.RData"))

print("Finished processing....")
q()
