suppressPackageStartupMessages({
  library(SCENIC)
  library(SCopeLoomR)
  library(scater)
})

args <- commandArgs(trailingOnly = TRUE)


#Directories
DIR_OUT <-"~/MBsnANA/HBana/SCENICana/SCENICOut/Class"
setwd(DIR_OUT)
Cls <- args[1]
Splt <- args[2]
## 1.--------
##load data
load(paste0(Cls,"/",Cls,"_",Splt,"_4loom.RData"))
SCE <- build_loom(paste0(Cls,"/",Cls,"_",Splt,".loom"),dgem=exprMat)
SCE <- add_cell_annotation(SCE, cellInfo)
close_loom(SCE)
print(paste0("finished creating the loom file...",Cls,"/",Cls,"_",Splt,".loom"))
print("Bye!")
dim(cellInfo)

q()