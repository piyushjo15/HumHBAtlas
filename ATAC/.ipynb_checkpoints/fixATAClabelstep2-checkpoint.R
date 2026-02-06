## THis script merged plotdata for ATAC from step2 RNA label transfer
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggrepel)
})

setwd("~/ATACana/Outs/")
## -----
load("plotdataHB_ATACv5.RData")
pdx <- plot.data.ATAC
cls <- unique(pdx$ClassX)
plot.data.ATAC <- c()
for(x in cls){
  load(paste0("Subset_",x,"/plotdata_RNAintstep2.RData"))
  pd <- data.frame(pd)
  plot.data.ATAC <- rbind(plot.data.ATAC,pd)
  rm(pd)
  
}
rm(x)
table(row.names(plot.data.ATAC)%in% row.names(pdx))

### also fix class now fow Neurobasts and Neurons
sub2cls <- read.delim("~/RNA/Classcluster.txt")
row.names(sub2cls) <- sub2cls$Cluster
head(sub2cls)
plot.data.ATAC$Class <- "ND"
clx <- unique(plot.data.ATAC$Cluster_step2)
for(x in clx){
  plot.data.ATAC[plot.data.ATAC$Cluster_step2==x,"Class"] <- sub2cls[x,"Class"]
}
table(plot.data.ATAC$Class)
cells <- row.names(plot.data.ATAC)
cells2 <- sample(cells,0.5*length(cells),replace = FALSE)
plot.data.ATAC$SplitTopic <- "B"
keep <- cells %in% cells2
plot.data.ATAC[keep,"SplitTopic"] <- "A"
##Fix cell IDs
cells <- row.names(plot.data.ATAC)
cells2 <- gsub("#","_",cells,fixed = TRUE)
cells2 <- gsub("-1","",cells2,fixed = TRUE)
plot.data.ATAC$CellID <- cells2
plot.data.ATAC$ATACID <- row.names(plot.data.ATAC)
save(plot.data.ATAC, file = "plotdataHB_ATACv6.RData")
