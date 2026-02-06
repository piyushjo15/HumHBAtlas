suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(batchelor)
  library(dplyr)
  library(tidyr)
  library(rlist)
  
})
### directories
DIROUT ="~/MBsnANA/HGGana/PFA/"
DIRS ="~/MBsnANA/HGGana/Scripts/RNA/"
setwd(DIROUT)
# 1. merging individual experiments ----
sams <- readLines(paste0(DIRS,"Gojo_PFA.txt"))
all.sce <- list()
all.dec <- list()
for(x in sams){
  all.sce[[x]] <- get(load(paste0(x,"scepro.RData")))
  all.dec[[x]] <- get(load(paste0(x,"dec.RData")))
  rm(list=ls(pattern="^P-"))

}

universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce <- lapply(all.sce, "[", i=universe,)
save(all.sce, file="comdata_GJPFA.RData")

universe2 <- Reduce(intersect, lapply(all.dec, rownames))
all.dec <- lapply(all.dec, "[", i=universe2,)
save(all.dec, file="comdatadec_GJPFA.RData")
#q()
## 2. plotdata ----

#load("comdata_GJPFA.RData")
#load("comdatadec_GJPFA.RData")

#plot data-----

sams <- names(all.sce)
plot.data <- c()
for(x in sams){
  del <- data.frame(colData(all.sce[[x]]))
  del2 <- reducedDim(all.sce[[x]], "UMAP")
  del$UMAP1 <- del2[,1]
  del$UMAP2 <- del2[,2]

  plot.data <- rbind(plot.data,del)
  rm(del,del2)
}

save(plot.data, file="plotdataGJPFA_121124.RData") ## fixed cell names and ind sample clustering, 310824
#rm(plot.data)
#q()
## 3. here I am trying to find 1op HVGs from each and combine them------
#load("plotdataGJPFA_101124.RData")
#not following combinevar
com.dec <- do.call(combineVar, all.dec)
com.hvg <- c()
#
sams <- names(all.dec)
for(x in sams){
  del <- getTopHVGs(all.dec[[x]])
  com.hvg <- c(com.hvg,del[1:1500])
}
com.hvg <- unique(com.hvg)
rm(all.dec)

load("~/MBnewana/HGNCsexchr.RData")
com.hvg <- com.hvg[!(com.hvg %in% SX)]

save(com.hvg, file="topHVGGJPFAnoRPMTSX_1500.RData")
#q()
# 
#load("topHVGGGJPFAnoRPMTSX_1500.RData")
## 4. cosine norm-------
sams <- unique(plot.data$Batch)
com.sce_cnnc <- c()
for(x in sams){
  del <- logcounts(all.sce[[x]])
  del <- cosineNorm(del[com.hvg,])
  com.sce_cnnc <- cbind(com.sce_cnnc,del)
  rm(del)
}
save(com.sce_cnnc, com.hvg, file="cosnormGJPFA.RData")
# q()
# rm(com.sce_cnnc)
## 5. combine counts----
com.sce <- do.call(noCorrect,list.append(all.sce, assay.type = "counts"))
assayNames(com.sce) <- "counts"
save(com.sce, file = "combinedcountsGJPFA.RData")
rm(com.sce)
## 6. MBN-------
norm.sce <- do.call(multiBatchNorm,
                    list.append(all.sce,norm.args=list(use_altexps=FALSE)))
rm(all.sce)

com.sce.rnnc <- do.call(noCorrect,norm.sce)
assayNames(com.sce.rnnc) <- "logcounts"

mdt <- logcounts(com.sce.rnnc)
save(mdt, file = "lgcountsGJPFA.RData")

print("Finished processing the script!!")
q()
