#This script is used to generate Arrow files from fragment files
## I have used my custome genome annotation derived from Gencode hg38 p13 r37
suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38) 
})

set.seed(456)
addArchRThreads(threads =1)  ##this only function properly with 1 thread
#define the Sample ID

sample = commandArgs(trailingOnly=TRUE)

##Dirs
DIR_ATACana <- "~/ATACana/"
DIR_FQ <- "~/snATACana/"
DIR_OUT <- "~/ATACana/Outs/Arrows/"

#Next, I load the reference genome
addArchRGenome("hg38")
load("extrafiles/GENCDH38p13r37_ann_4_archr.RData") ##custom gene annotation 

#define the input file
inputFiles <- paste0(DIR_FQ,sample,"/outs/fragments.tsv.gz")

#generate ArrowFiles
setwd(DIR_OUT)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sample,
  geneAnnotation=new_gen_ann,
  minTSS = 3,  
  minFrags = 3000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = "001_barcode_qc",
  promoterRegion = c(2000,100),
  force = T
)
print("Generated Arrow File ")


q()
