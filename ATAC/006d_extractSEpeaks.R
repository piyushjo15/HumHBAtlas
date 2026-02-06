suppressPackageStartupMessages({
  library(tidyverse)
})


setwd("~//ATACana/Outs/CRI/rose")
cls <- commandArgs(trailingOnly = TRUE)

## -------
## maybe it is better to obtain a list per class and then select top 100 from that
## For each of the SE obtian a consistuent CREs and make a list of those CREs

sebed <- read.delim(paste0(cls,"_SE.bed"),header = FALSE)
sebed <- sebed[order(sebed$V4),]
sebed$size <- sebed$V3-sebed$V2
sebed$Signal <- sebed$V4/sebed$size

## remove SEs that are less than 5000 kb in size
sebed <- sebed[sebed$size>4999,]
sebed$Rank <- rank(sebed$Signal,ties.method = "first")
sebed <- sebed[order(sebed$Signal,decreasing = TRUE),]
## take at max 150
nses <- dim(sebed)[1]
if(nses>150){
  sebed <- sebed[1:150,]
}
#head(sebed)
#tail(sebed)

sebed$CRE <- paste0(sebed$V1,":",sebed$V2,"-",sebed$V3)

bed <- read.delim("~/ATACana/Outs/PeakCom/AllPeaks_robust.fil.sorted.bed",header = FALSE)

#head(bed)

se_cre_list <- list()
cres <- sebed$CRE
for(i in 1:length(cres)){
  chr <- sebed[i,"V1"]
  a <- sebed[i,"V2"]
  b <- sebed[i,"V3"]
  bedx <- bed[bed$V1==chr,]
  bedx <- bedx[bedx$V2>=a & bedx$V3<b,"V4"]
  se_cre_list[[cres[i]]] <- bedx
  rm(chr,a,b,bedx)
  
}
rm(i)

save(se_cre_list, file = paste0(cls,"_SE_cre.RData"))
q()

## part 2 -----------
## Finding overlap among SE CREs across classes
Cls <- readLines("~/RNA/RNA_Class.txt")
all_se_cres <- list()
metase <- c()
for(x in Cls){
  load(paste0(x,"_SE_cre.RData"))
  secres <- names(se_cre_list)
  del <- data.frame(SE=secres,ID=paste0("SE_",x,"_",seq(length(secres))),Class=x)
  metase <- rbind(metase,del)
  names(se_cre_list) <- del$ID
  all_se_cres <- do.call(c, list(all_se_cres, se_cre_list))
  rm(se_cre_list,secres,del)
}
rm(x)

all_ses <- row.names(metase) <- metase$ID
metase$Len <- lengths(all_se_cres)
save(all_se_cres,metase, file = "ALLSElist.RData")
