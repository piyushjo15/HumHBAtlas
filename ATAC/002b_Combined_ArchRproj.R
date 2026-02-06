## This script is used for a combined project
## This is for copying original unaltered arrows files into a new folder
## for integrated analysis of all samples together.
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyr)
  library(dplyr)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})


addArchRThreads(threads =1)
#Next, I load the reference genome
addArchRGenome("hg38")

DIR="/home/p541i/MBsnANA/HBana/ATACana/Outs/"
DIRO="/home/p541i/MBsnANA/HBana/ATACana/Outs/comATAC/"
setwd(DIR)
## cannot identify which SA283.arrow file I uses
## most probably the one located in SA283 based on output of 001c_extraprocessing.R
#Define Arrow files
Arrow_Path <- paste0(DIR,"original_arrows/")
ArrowFiles <- list.files(Arrow_Path)
ArrowFiles <- paste0(Arrow_Path,ArrowFiles)

proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = DIRO,
    copyArrows = TRUE
)

###subset to filtered cells
load("extrafiles/plotdataHB_ATAC.RData")
cn <- row.names(plot.data.ATAC)
table(cn %in% getCellNames(proj))
proj <- proj[cn,]
#saveArchRproj in between
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = DIRO)

q()
