## This scripts converts the RNA RData into python anData, that will be used as 
## SCENIC+ input. 
## This is the only Group 3/4 tumor sample that is not truly multi-omic 
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/path/to/python3.7")
use_python("~/path/to/python3.7")

suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(SingleCellExperiment)
  library(BSgenome.Hsapiens.UCSC.hg38)
})


#load RNA SCE
Sample <- "SA194"
RNA_Sample <- "SN407"

DIR_RNA <- "RNA/"


##loading log normalized RNA data, no imputed cells in this sample
seRNA <- get(load(paste0(DIR_RNA,Sample,"/",Sample,"_SCE.RData"))) 

## combined tumor single-cell HVG list
com.hvg <- readLines("Files/topHVGG34new10xnoRPMTSX.txt") #com.hvg

#convert RNA SCE to Scanpy AnnData
############ add missing information and filter for ATAC cells ###############

#remove normal cells from sceRNA
remove_cells <- readLines(paste0(DIR_RNA,"normal_cells/",Sample,"_normal_cells.txt"))

print("Number of normal cells")
length(normal_cells)

print("Number of cells before removing Normal cells")
length(colnames(seRNA))

seRNA <- seRNA[,colnames(seRNA) %ni% normal_cells]

print("Number of cells after removing Normal cells")
length(colnames(seRNA))


######## subset RNA genes to contain only HVGs ############
#filter for HVG
#subst the SCE to contain only the top HVG and the cells also present in the ATAC project
seRNA <- seRNA[rownames(seRNA) %in% com.hvg,]

#save seRNA used for SCENIC+
print("Saving seRNA object used for SCENIC+")
save(seRNA, file = paste0(DIR_RNA,Sample,"/",Sample,"_seRNA_SCENIC.RData"))


############ cp the ind_cluster as a new column naming predictedGroup (to have the same name as in ATAC data)

seRNA$predictedGroup <- seRNA$ind_cluster


############# step-by-step to Anndata ###############

#save the matrix
a <- seRNA@assays@data$logcounts
colnames(a) <- colnames(seRNA)

#save the barcodes
barcodes<-data.frame(colnames(seRNA))
colnames(barcodes)<-'Barcode'

#Save the gene names
genes<-data.frame(rownames(seRNA))
colnames(genes)<-'Gene'

#save metadata
cellMeta<-as.data.frame(seRNA@colData)

######### Save the data to python ############
repl_python()

import scanpy as sc
import pandas as pd
from scipy import io
import os

Sample = r.Sample
projDir = os.path.join("SCENIC/output/RNA/",Sample)

#log count Matrix 
counts = r.a
from scipy import sparse
counts=sparse.csr_matrix(counts)

#load metadata
barcodes=r.barcodes
genes=r.genes

adata=sc.AnnData(counts.T)
adata.raw = adata
adata.obs_names=barcodes['Barcode'].values
adata.var_names=genes['Gene'].values

#import metadata
cellMeta=r.cellMeta
adata.obs=cellMeta


#save adata object
adata.write(os.path.join(projDir, 'adata.h5ad'), compression='gzip')

exit
q()
print(paste0("Saved adata object for sample: ", Sample))
