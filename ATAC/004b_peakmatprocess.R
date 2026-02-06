## getting peaks called on cluster_subset
suppressPackageStartupMessages({
  library(ArchR)
   library(tidyverse)
})

set.seed(456)
addArchRThreads(threads = 6)
cls <- commandArgs(trailingOnly = TRUE) ## this has to be done per class else too much memory
#setwd
#Next, I load the reference genome
addArchRGenome("hg38")
DIR <- "~/ATACana/Outs/"
setwd(DIR)

##get peak matrix per cluster-----
proj <- loadArchRProject(path = "comATACx/")
df <- data.frame(getCellColData(proj))
cells <- row.names(df)[df$predictedGroup==cls]
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
#save(PeakMatrix, file = paste0("comATACx/",cls,"_PeakMatrix.RData"))
region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)
PM <- assay(PeakMatrix)
region_names$Freq  <- rowSums(PM > 0)/length(cells)
colnames(region_names)[7] <- paste0(cls,"_Freq")
save(region_names, file=paste0("PeakCom/",cls,"_proppeask.RData"))

#q()
### filtering Peaks ######
cls <- readLines("~/RNA/Class.txt")
df <- get(load(paste0("PeakCom/",cls[1],"_proppeask.RData")))
head(df)
for(i in 2:length(cls)){
  load(paste0("PeakCom/",cls[i],"_proppeask.RData"))
  df <- cbind(df, region_names[7])
}
head(df)
df2 <- df[,-c(1:6)]
head(df2)
PS <- data.frame(PeakSet)
row.names(PS) <- paste0(PS$seqnames,":",PS$start,"-",PS$end)
table(row.names(PS)==row.names(df2))
PS$max_freq <- apply(df2, 1, max)
sum(PS$max_freq >= 0.01)
sum(PS$max_freq >= 0.03)
sum(PS$max_freq >= 0.05)
PS$robust <- PS$max_freq >= 0.03 
save(PS,  file = "PeakCom/Peaks_robust.RData")
bed <- PS[PS$robust,c("seqnames",  "start",    "end")]
head(bed)
write.table(bed, file="PeakCom/AllPeaks_robust.bed", row.names = FALSE,col.names = FALSE, quote = FALSE, sep="\t")
##these peaks are sorted later
sort -k1,1 -k2,2n AllPeaks_robust.bed > AllPeaks_robust.sorted.bed
bedtools subtract -A -a AllPeaks_robust.sorted.bed -b ~/extrafiles/hg38-blacklist.v2_noncanchr.bed > AllPeaks_robust.fil.sorted.bed

