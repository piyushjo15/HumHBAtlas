library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")


suppressPackageStartupMessages({
  library(Matrix)
})
##
setwd("~/braun")

# ## part 1----------
# ## First extract the metadata and gene ids to check how to subset and so on
# 
# repl_python()
# import os
# import scanpy as sc
# import numpy as np
# import pandas as pd
# from scipy.sparse import csc_matrix
# 
# adata = sc.read_h5ad("Braun.h5ad")
# print("here is the data description...")
# adata
# print('here are the keys to layers...')
# adata.layers.keys()
# df=adata.obs
# genes=adata.var
# exit
# pd <- py$df
# genes <- py$genes
# save(pd,genes, file = "Braun_meta.RData")
# q()
# 
# ## part 2---------
# ## now check genes and how to subset data
#load("Braun_meta.RData")
#pd$assay <- as.character(pd$assay)
#pd <- pd[pd$assay=="10x 3' v3",]
## get cells per batch so I can logtransform per batch, that seems to be the easiest way for now
#pd$sample_id <- as.character(pd$sample_id)
#btc <- unique(pd$sample_id)
#length(btc)
#head(btc)
#class(btc)
#for(x in btc){
#    cells <- row.names(pd)[pd$sample_id==x]
#    write(cells, file=paste0("cellid_",x,".txt"))
#}
#write(btc, file="sampleids_braunv3.txt")
#
##fix cell meta-----
# ##remove 10xv2
# pd$assay <- as.character(pd$assay)
# 
# pd <- pd[pd$assay=="10x 3' v3",]
# 
# pd$CellClass <- as.character(pd$CellClass)
# pd$sample_id <- as.character(pd$sample_id)
# pd$cluster_id <- as.numeric(pd$cluster_id)
# pd$Region <- as.character(pd$Region)
# 
# ##now merge erythrocyte with immune
# pd[pd$CellClass=="Erythrocyte","CellClass"] <- "Immune"
# 
# ##now merge Neural crest and vascular with fibroblast
# pd[pd$CellClass=="Neural crest","CellClass"] <- "Fibroblast"
# pd[pd$CellClass=="Vascular","CellClass"] <- "Fibroblast"
# 
# ##remove placodes
# pd <- pd[pd$CellClass!="Placodes",]
# save(pd, file = "Braun_metav3.RData")
# pd1 <- pd
# ## subsetting to max 200 cells per cluster
# 
# pd$cluster_id <- paste0("Cl",pd$cluster_id)
# 
# cells <- c()
# 
# cls <- unique(pd$cluster_id)
# for(x in cls){
#   del <- row.names(pd)[pd$cluster_id==x]
#   if(length(del)>200){
#     del <- sample(del,200)
#   }
#   cells <- c(cells,del)
# 
# 
# }
# pd <- pd[row.names(pd) %in% cells,]
# 
# save(pd, file = "Braunmeta_subsettedbycluster.RData")
# 
# ## subsetting to max 2000 cells per supercluster
# pd <- pd1
# 
# cells <- c()
# 
# cls <- unique(pd$CellClass)
# for(x in cls){
#   del <- row.names(pd)[pd$CellClass==x]
#   if(length(del)>2000){
#     del <- sample(del,2000)
#   }
#   cells <- c(cells,del)
# 
# 
# }
# pd <- pd[row.names(pd) %in% cells,]
# 
# save(pd, file = "Braunmeta_subbyclass.RData")
# 
# # ## fix genes meta --------
# ens <- read.delim("~/extrafiles/gencdH38p13r37CR_genes.txt",header = FALSE)
# ens <- ens %>% separate(V1,c("A","B"))
# keep <- duplicated(ens$A)
# ens <- ens[!keep,]
# row.names(ens) <- ens$A
# ex <- row.names(genes)
# table(ex %in% row.names(ens))
# ex <- ex[ex %in% row.names(ens)]
# genes <- genes[ex,]
# genes$Gene2 <- ens[ex,"V2"]
# save(genes, file = "filteredgenes_braun.RData")
# write(row.names(genes), file = "filter_ensg_braun.txt")

## part 3----------
## obtain logcounts for selected cells
suppressPackageStartupMessages({
  library(scran)
  library(scater)
})

sams <- readLines("~/braun/sampleids_braunv3.txt")
load("filteredgenes_braun.RData")
#load("Braunmeta_subbyclass.RData")
load("Braunmeta_subsettedbycluster.RData")

## to make sure genes are correctly ordered
genes <- genes$Gene2
pd$sample_id <- as.character(pd$sample_id)
head(pd)
head(genes)

mdt <- c()


for(id in sams){
    load(paste0("lgmdt_",id,".RData"))
    cells <- row.names(pd)[pd$sample_id==id]
    mat_del <- mat_del[genes,cells]
    mdt <- cbind(mdt,mat_del)
    rm(mat_del,df,cells)
}
table(is.na(mdt))

#save(pd,mdt, file = "Braun_proc_subbyclass.RData")
save(pd,mdt, file = "Braun_proc_subbycluster.RData")


q()


####







##