## This script aguments the region in a bed, resulting in 5 overlapping regions per region
## extended to 100 bp in both direction by a stride of 50 bp
suppressPackageStartupMessages({
  library(tidyverse)
})

DIRO <- "~/DeepHB/Outs/"

setwd(DIRO)
args <- commandArgs(trailingOnly = TRUE)
#

##function
augmentation_shifts <- c((-200),(-150),(-100), (-50), 0, 50, 100,150,200)

augment_peak <- function(df, shifts=augmentation_shifts) {
  ## Requires a single peak
  if (nrow(df) != 1) {
    stop("This function expects a single peak as an input")
  }
  
  new_df <- data.frame(
    chr=df$chr,
    start=df$start + shifts,
    end=df$end + shifts,
    peak=df$peak,
    aug_peak=paste0(df$peak, "_", shifts),
    aug_shift=shifts
  )
  return(new_df)
}

## for

#inp_bed <- read.delim("PedCanGWAS/orig_bed/all.bed",header = FALSE)
inp_bed <- read.delim(paste0("HGG/H3K27acDEP/",args,"_marker.bed"),header = FALSE)
inp_bed$V4 <- paste0(inp_bed$V1,":",inp_bed$V2,"-",inp_bed$V3)
#keep <- is.na(inp_bed$V4)
#k <- grep("TRUE",keep)
#inp_bed[keep,"V4"] <- paste0("rsND",seq(length(k)))
#inp_bed$V5 <- inp_bed$V4

#colnames(inp_bed) <- c("chr","start","end","peak","rsid")
colnames(inp_bed) <- c("chr","start","end","peak")

head(inp_bed)

out_bed <- Reduce(bind_rows,lapply(unique(inp_bed$chr), function(x) {
  chr_df <- filter(inp_bed, chr==x)
  chr_df_aug <- Reduce(bind_rows, lapply(1:nrow(chr_df), function(i) {
    augment_peak(chr_df[i,])
  }))
  return(chr_df_aug)
})) %>%
  left_join(select(inp_bed, -chr, -start, -end))


head(out_bed)
#write.table(out_bed, file = "PedCanGWAS/all.aug.bed",
#            sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(out_bed, file = paste0("HGG/H3K27acDEP_aug/",args,"_marker.bed"),
            sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
q()
###
