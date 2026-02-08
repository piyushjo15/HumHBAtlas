## This script aguments the region in a bed, resulting in 5 overlapping regions per region
## extended to 100 bp in both direction by a stride of 50 bp

suppressPackageStartupMessages({
  library(tidyverse)
})

DIRO <- "~/DeepHB/"

setwd(DIRO)
args <- commandArgs(trailingOnly = TRUE)
print(paste0("Augmenting peaks for topic MP:",args))
#

##function
augmentation_shifts <- c((-100), (-50), 0, 50, 100)

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
## These bed files are obtained from meta-regulatory program (MP) analysis of ATAC data
## each MP is represented by a set of 2500 peak regions
## a radon set of 2500 CRE was also used for backgroun MP00.
inp_bed <- read.delim(paste0("TopicMP_CREs/",args,"_CREs.bed"),header = FALSE) 

inp_bed <- inp_bed[,c(1:4)]

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


#head(out_bed)
write.table(out_bed[, c("chr", "start", "end", "aug_peak")],
            paste0("TopicMP_CREs_aug/",args, "_CREs_aug.bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)

q()

