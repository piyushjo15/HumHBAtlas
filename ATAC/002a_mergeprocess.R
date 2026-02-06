## Here I am merging the metadata per library and also merging gene-score matrices
## to obtain a combined gene-score matrix (GSM) for ATAC data
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(rlist)
})
DIR_ATAC <- '~/ATACana/Outs/'
setwd(DIR_ATAC)

## merging individual experiments GSM----
libs <- readLines("HBlibsatac.txt")

plot.data.ATAC <- c()
for(x in libs){
 load("projs/",x,"/",x,"_barcode_stats.RData")
  plot.data.ATAC <- rbind(plot.data.ATAC,plot.data)
 rm(plot.data)
}
rm(x)
plot.data.ATAC$ATAC_barcodes <- row.names(plot.data.ATAC)

save(plot.data.ATAC, file="plotdataHB_ATAC.RData")

## merging individual experiments GSM----
 all.sce.atac <- list()
##just lopping and rbinding is not working, need to use scran/scater functions
 
for(x in libs){
 load(paste0(DIR_ATAC,"projs/",x,"/",x,"_GSM.RData"))
  sce <- SingleCellExperiment(list(counts=mdt))
  cl2 <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = cl2, min.mean = 0.1 )
  sce <- logNormCounts(sce)
  all.sce.atac[[x]] <- sce
  rm(mdt,sce)
}
rm(x)
save(all.sce.atac, file = paste0(DIR_ATAC,"combined_sce_ATAC.RData"))
#combine counts-------
#load("combined_sce_ATAC.RData")
com.sce <- do.call(noCorrect,list.append(all.sce.atac, assay.type = "counts"))
mdt.counts <- assay(com.sce)
save(mdt.counts, file = paste0(DIR_ATAC,"merged_ATACGSM.RData"))
rm(com.sce,mdt.counts)
##logcounts
#combine counts
com.sce.nrnc <- do.call(noCorrect,list.append(all.sce.atac, assay.type = "logcounts"))
mdt.lgcounts <- assay(com.sce.nrnc)
save(com.sce, file = "merged_ATACGSM.RData")

q()
