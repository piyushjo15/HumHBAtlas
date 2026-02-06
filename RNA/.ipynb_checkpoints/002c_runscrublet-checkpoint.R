## This script takes SoupX corrected data to predict doublets
## This requires python
library(reticulate)
##python
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")
#argument
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(Matrix)
})

args = commandArgs(trailingOnly=TRUE)
load(paste0(args[1],"_postSoupXDX.RData"))

#scrublet----
mdt <- counts(sce)
mdt <- t(mdt)
# python code
repl_python()
import scrublet as scr
scrub = scr.Scrublet(r.mdt)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
exit
sce$SCRscore <- py$doublet_scores
keep <- py$predicted_doublets
SCRcall <- rep("Singlet", length(keep))
SCRcall[keep] <- "Doublet"
sce$SCRcall<- SCRcall

sce <- subset(sce, ,SCRcall=="Singlet")
save(sce, mdfr2,mdt.pd, file = paste0(args,"_postSCR.RData"))
q()
