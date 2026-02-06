## Processed data from each library was merged together to obtain
## a combined metadata and gene-expression data
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(rlist)
})
libs <- read.delim("all_libs.txt",header = FALSE)
libs <- libs$V1
## merging individual experiments ----
all.sce <- list()
all.dec <- list()
# # plotdata
cl <- c("nUMIs","nGenes","pct.mt","MALAT1","in_ex_frac","decontX_contamination",
        "SCRscore","sizeFactor","total_counts","S.score","G2M.score","CC.score",
        "CC.per","ind_cluster","Batch","Stage", "Age","Sex",
        "Sample","Region","iUMAP1","iUMAP2")

plot.data <- c()
for(x in libs){
  sce <- get(load(paste0(x,"scepro.RData")))
  del <- data.frame(colData(sce))
  del2 <- reducedDim(sce, "UMAP")
  del$iUMAP1 <- del2[,1]
  del$iUMAP2 <- del2[,2]
  all.sce[[x]] <- sce
  plot.data <- rbind(plot.data,del[,cl])
  rm(del,del2,sce)

  all.dec[[x]] <- get(load(paste0("files/",x,"dec.RData")))
  rm(list=ls(pattern="^SN"))

}
rm(x)
save(plot.data, file="plotdataHB.RData")

universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce <- lapply(all.sce, "[", i=universe,)
save(all.sce, file="comdata.RData")

universe2 <- Reduce(intersect, lapply(all.dec, rownames))
all.dec <- lapply(all.dec, "[", i=universe2,)
save(all.dec, file="comdatadec.RData")



#here I am trying to find 1op HVGs from each and combine them------
#load("comdatadec.RData")

##hvg using combine var
com.dec <- do.call(combineVar,all.dec)
topHVG <- getTopHVGs(com.dec)

load("extrafiles/HGNCsexchr.RData")
topHVG <- topHVG[!(topHVG %in% SX)]
rm(all.dec)

save(topHVG, file="topHVGcomhvgnoRiboMTSX.RData")
# #q()
##cosine norm----
#load("topHVGcomhvgnoRiboMTSX.RData")
com.sce_cnnc <- c()
for(x in libs){
  sce <- get(load(paste0(x,"scepro.RData")))
  rm(list=ls(pattern="^SN"))

  #sce <- all.sce[[x]]
  del <- logcounts(sce)
  del <- cosineNorm(del[topHVG[1:4000],]) 
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del,sce)
}
rm(x)
save(com.sce_cnnc, topHVG, file="cosnormHB.RData")
dim(com.sce_cnnc)
rm(com.sce_cnnc)
q()
# ##------
load("comdata.RData")
#combine counts
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsHB.RData")
#q()
rm(com.sce)


norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce.rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce.rnnc) <- "logcounts"
mdt <- logcounts(com.sce.rnnc)
save(mdt, file = "lgcountsHB.RData")

q()