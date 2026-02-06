## Getting Index file per sample, each file has a column for cell id and a column
## for annotation
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
})

DIR <- "~/ATACana/Outs/"
DIRSE <- "~/ATACana/Outs/CRI/"
DIRSEB <- "~//ATACana/Outs/CRI/BGD/"

setwd(DIR)

## for each Class I will extract Cells for each sample
load("plotdataHB_ATACv6.RData")
plot.data.ATAC$Class_split <- paste0(plot.data.ATAC$Class,"_",plot.data.ATAC$SplitTopic)
#subset further, looks like 3k cells should be good
cells <- c()
cls <- unique(plot.data.ATAC$Cluster_step2)
for(x in cls){
  cells_sel <- row.names(plot.data.ATAC)[plot.data.ATAC$Cluster_step2==x]
  if(length(cells_sel)>3000){
    cells_sel <- sample(cells_sel,3000)
    cells <- c(cells,cells_sel)
  } else {
    cells <- c(cells,cells_sel)
  }
  rm(cells_sel)
}
rm(x)
plot.data.ATAC <- plot.data.ATAC[cells,]

## classes that are affected by subsetting
## GC, GCUBC, Glioblasts, Immune, MixedGlia, MOG, Neuroblasts_II,
## Neuroblasts_III, Neurons_I, Neurons_II, PreO, Pro and PC
## Only MES_VAS, CBINT and Neuroblasts_I are not downsampled

df <- plot.data.ATAC[,c('Sample',"Class_split")]
df$Index <- row.names(df)
df <- df %>% separate(Index,c("A","B"), sep="#")
head(df)
clss <- unique(df$Class_split)
for(x in clss){
  del <- df[df$Class_split==x,]
  sam_cnt <- data.frame(table(del$Sample))
  ##sample with more than 49 cells
  sams <- as.character(sam_cnt$Var1)[sam_cnt$Freq>49]
  del <- del[del$Sample%in%sams,]
  paths=paste0(DIRSE,x)
  if(!dir.exists(paths)){
    dir.create(paths, showWarnings = TRUE, recursive = FALSE)
  }
  write(sams, file = paste0(DIRSE,x,"/Samples.txt"))
  #head(del)
  for(y in sams){
    del2 <- del[del$Sample==y,"B"]
    write(del2, file = paste0(DIRSE,x,"/",y,"_Index.csv"))
    rm(del2)
  }
  rm(del,sams,y)
}

rm(x)
## background -------
## I am randomly selecting 10k cells to generate a control bam
## for superenhancer analysis
load("plotdataHB_ATACv6.RData")
rn <- row.names(plot.data.ATAC)
rnx <- sample(rn,10000,replace = FALSE)
pd <- plot.data.ATAC[rnx,]
pd$Cell_type <- "Control"
pd$Index <- row.names(pd)
pd$ATAC_barcodes <- NULL
pd <- pd %>% separate(Index,c("A","B"), sep="#")
head(pd)
sam_cnt <- data.frame(table(pd$Sample))
##sample with more than 49 cells
sams <- as.character(sam_cnt$Var1)[sam_cnt$Freq>49]
pd <- pd[pd$Sample%in%sams,]

if(!dir.exists(DIRSEB)){
  dir.create(DIRSEB, showWarnings = TRUE, recursive = FALSE)
}
write(sams, file = paste0(DIRSEB,"Samples.txt"))
#head(del)
for(y in sams){
  del2 <- pd[pd$Sample==y,"B"]
  write(del2, file = paste0(DIRSEB,"/",y,"_Index.csv"))
  rm(del2)
}
rm(del,sams,y)


q()