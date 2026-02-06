### Subsetting ATAC clusters by groups to improve clustering
## Smaller and poor quality clusters obtained after this clustering were remvoed
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(foreach)
  library(tidyverse)
})

args = commandArgs(trailingOnly=TRUE) ## Lv1 grouping

set.seed(456)
addArchRThreads(threads = 8)

#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/ATACana/Outs/"
setwd(DIR)
Sample="comATAC/"


################ ATAC: load the data ########################
proj <- loadArchRProject(path = Sample)
###sub-setting by  super clusters
plot.data <- data.frame(getCellColData(proj))
## the 'Clusters_int' called in 002c_CombinedATAC_LSI_UMAP.R were merged into 11 groups
## and these groups will be subsetted to improve clustering resolution
cl_lv <- read.delim("Clusters_int.txt",row.names=1) 
keep <- cl_lv$Supercluster==args
sel <- row.names(cl_lv)[keep]
keep <- plot.data$Clusters_int %in% sel
table(keep)
cells <- row.names(plot.data)[keep]
rm(plot.data)
print("Subsetting ArchRproject by supercluster:")
subsetArchRProject(proj,
                   cells=cells,
                   outputDirectory = paste0("Subset_",args))
rm(proj)

### subclustering and increasing clustering rsolution ------
proj <- loadArchRProject(paste0("Subset_",args))

print("Adding clusters:")
proj <- addClusters(input = proj,
                    name = "Clusters_subset",
                    reducedDims = "IterativeLSI_int",
                    method = "Seurat",
                    force = T,
                    maxClusters = 50, #increased from 25
                    resolution=2, ##increased resolution
                    corCutOff = 0.75,
                    scaleDims = FALSE,
                    seed = 1)

print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = paste0("Subset_",args))

print("Saving metadata:")
plot.data <- data.frame(getCellColData(proj))
save(plot.data, file = paste0("Subset_",args,"/metadata_SC.RData"))
### i merged individual metadata into plotdataHB_ATACv3.RData
q()
