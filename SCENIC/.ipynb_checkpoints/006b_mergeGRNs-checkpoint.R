##This script merges GRNs from all the classes

suppressPackageStartupMessages({
  library(tidyverse)
})
DIR = "~/SCENICOut/Class/"
#read metadata
setwd(DIR)
Cls <- readLines("NA_Class.txt")
##This script merges GRNs from all the classes
### collecting all the GRNs identified in all classes.
class_grn_list <- list()
all_tfs <- c()
for(x in Cls){
  load(paste0(x,"/SCENIC/",x,"_gmt.RData"))
  class_grn_list[[x]] <- regs
  tfdel <- names(regs)
  all_tfs <- rbind(all_tfs,data.frame(TF=tfdel,Class=rep(x,length(tfdel))))
  rm(regs,tfdel)
  
  
}
tfs <- unique(all_tfs$TF)
all_grns <- list()
for( x in tfs){
  del <- all_tfs[all_tfs$TF==x,"Class"]
  for(y in del){
    del_grn <- class_grn_list[[y]]
    all_grns[[paste0(x,"_",y)]] <- del_grn[[x]]
    rm(del_grn)
  }
  rm(del)
}
lengths(all_grns)
save(all_grns,all_tfs, file = "AllGRNs_v2.RData")
q()
