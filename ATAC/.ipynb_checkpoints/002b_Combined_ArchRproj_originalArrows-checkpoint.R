## This script is used for a combining arrow files for all the ATAC data into 
## single ArchR  project for integrated analysis of all samples together.
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

addArchRThreads(threads =1)
#Next, I load the reference genome
addArchRGenome("hg38")

DIR="~/ATACana/Outs/"
setwd(DIR)
#Define Arrow files
Arrow_Path <- paste0(DIR,"Arrows/") ##arrow files were also copied for backup
ArrowFiles <- list.files(Arrow_Path)
ArrowFiles <- paste0(Arrow_Path,ArrowFiles)

proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = "comATAC",
    copyArrows = TRUE
)

###subset to filtered cells
load("plotdataHB_ATAC.RData")
cn <- row.names(plot.data.ATAC)
table(cn %in% getCellNames(proj))
proj <- proj[cn,]
#saveArchRproj in between
print("Saving ArchR project:")
saveArchRProject(proj, outputDirectory = "comATAC")

q()
