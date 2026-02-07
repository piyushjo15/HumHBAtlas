##This script fixes GRN per class and split, and then merges GRNs 
## obtained per class from both splits

suppressPackageStartupMessages({
  library(tidyverse)
})
#get a TF to region matrix
Cls <- commandArgs(trailingOnly = TRUE)
DIR = "~/SCENICOut/Class/"
#read metadata
setwd(paste0(DIR, Cls))

## part 1 ------------------------
## First part is to extract predicted GRNs from each split run and fix them
## Fixing is this. Previously I had seen some TF to region and region to target
## links but not TF to region to target links. This was becuase the targets
## were dropped due to some ranking issue. I fixed this.
Splt <- "A" # A or B
## v2 is using importance_x_rho for tf to gene rankings that saved OLIG2 eGRNs in PreO
## The out of v2 input has "_v2" as suffix
df <- read.csv(paste0("SCENIC/",Cls,"_",Splt,"/GRNs/",Cls,"_",Splt,"_eRegulon_metadata.csv"), sep = ";")
df2 <- df[df$is_extended!="True",] ## In G34MB, we used is_extended==True for some reasons
gg <- as.character(df2$Gene_signature_name)
ggx <- grep("+_+",gg, fixed = TRUE) ## only selecting positive regulation
## "+_-" are also positive TF-target relations with enhancer closed so I am not
## using them

##extracting enhancer identified in this sample
df3 <- df2[ggx,]
tfs <- unique(df3$TF)
rm(df,df2)
##unique enhancer regions
reg_adj <- data.frame(Region=df3$Region)
reg_adj <- reg_adj %>% separate(Region,c("A","B","C"))

keep <- duplicated(df3$Region)
table(keep)
reg_adj <- reg_adj[!keep,]

head(reg_adj)
##this bed needs to be sorted
write.table(reg_adj, file=paste0("SCENIC/",Cls,"_",Splt,"/GRNs/",Cls,"_",Splt,"_enhancers.bed"),
            sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

###extracting regions associated with a TF for positive regulations
TF_region <- list()
for (x in tfs){
  del <- df3[df3$TF==x,c("TF","Region")]
  del <- del[!duplicated(del$Region),] ##remove duplicates
  TF_region[[x]] <- del
  rm(del)
}
rm(x)
#get a region to gene matrix
regs <- unique(df3$Region)
region_tar <- list()
for (x in regs){
  del <- df3[df3$Region==x,c("Region","Gene", "R2G_importance", "R2G_rho", "R2G_importance_x_rho")]
  XX <- paste0(del$Region,"_",del$Gene)
  del <- del[!duplicated(XX),]#remove duplicates
  region_tar[[x]] <- del
  rm(del,XX)
}
rm(x)

##now I want three things, TF, regions and gene links
GRNS <- c()
for(x in tfs){
  regions <- as.character(TF_region[[x]]$Region)
  targets <- c()
  for (y in regions){
    del <- region_tar[[y]]
    targets <- rbind(targets,del)
    rm(del)
  }
  targets$TF <- x
  GRNS <- rbind(GRNS,targets)
  rm(y,targets,regions)
}
rm(x)

table(duplicated(paste0(GRNS$Region,"_",GRNS$Gene,"_",GRNS$TF))) ##checking duplicated connections, ALL false

row.names(GRNS) <- seq(dim(GRNS)[1])
head(GRNS)
#get a region to gene matrix
DIRx <- paste0("~/MBsnANA/HBana/SCENICana/SCENICOut/Class/",Cls)
TF2G <- read.delim(paste0(DIRx,"/cor",Cls,"_fil.tsv")) 

colnames(TF2G)[c(3,5)] <- c("TF2G_importance","TF2G_rho")
keep <- is.na(TF2G$TF2G_rho)
TF2G <- TF2G[!keep,]
##selecting only pos corr
TF2G <- TF2G[TF2G$TF2G_rho>=0.05,] ##adding cut-off rho
GRNS2 <- c()
## just addition check
tfs <- tfs[tfs %in% TF2G$TF]

for(x in tfs){
  TF2G_sel <- TF2G[TF2G$TF==x,]
  ##add TF-TF linksif not present
  if(!(x %in% TF2G_sel$target )){
    imp <- max(TF2G_sel$TF2G_importance)+ 1e-05
    df_del <- data.frame(x,x,imp,1,1)
    colnames(df_del) <- colnames(TF2G_sel)
    TF2G_sel <- rbind(df_del,TF2G_sel)
    rm(df_del)
  }
  tars1 <- unique(TF2G_sel$target)
  GRNS_sel <- GRNS[GRNS$TF==x,]
  tars2 <- unique(GRNS_sel$Gene)
  targets <- intersect(tars1,tars2)
  ##subset to shared targets
  GRNS_sel <- GRNS_sel[GRNS_sel$Gene %in% targets,]
  GRNS_sel$TF2G_importance <- GRNS_sel$TF2G_rho <- NA

  for (y in targets){
    GRNS_sel[GRNS_sel$Gene==y,"TF2G_importance"] <- round(TF2G_sel[TF2G_sel$target==y,"TF2G_importance"],4)
    GRNS_sel[GRNS_sel$Gene==y,"TF2G_rho"] <-round(TF2G_sel[TF2G_sel$target==y,"TF2G_rho"],4)
  }
  GRNS2 <- rbind(GRNS2,GRNS_sel)
  rm(y,targets,GRNS_sel,TF2G_sel,tars1, tars2 )
}
rm(x)
GRNS2 <- GRNS2[,c(6,1:5,7,8)]
row.names(GRNS2) <- seq(dim(GRNS2)[1])
tail(GRNS2)
write.table(GRNS2,paste0("SCENIC/",Cls,"_",Splt,"/GRNs/",Cls,"_",Splt,"_eRegulons_fix.txt"),
            quote = FALSE, row.names = FALSE, sep ="\t")

##getting TF-GRN
def <- GRNS2[,c("TF","Gene")]
tfs <- unique(def$TF)
regs <- list()
for(x in tfs){
  defx <-def[def$TF==x,]
  regs[[x]] <- unique(c(defx$Gene,x))
  rm(defx)
}
rm(x)
save(regs, file = paste0("SCENIC/",Cls,"_",Splt,"/GRNs/",Cls,"_",Splt,"_gmt.RData"))
#make a bed file
GRNS2 <- GRNS2 %>% separate(Region,c("chr","start","end"))
head(GRNS2)
GRNS2$ID <- paste0(GRNS2$TF,"_",GRNS2$start,"_",GRNS2$Gene)
bed <- GRNS2[,c("chr","start","end","ID")]
head(bed)
write.table(bed,paste0("SCENIC/",Cls,"_",Splt,"/GRNs/",Cls,"_",Splt,"_enhancersv2.bed"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep ="\t")

##extracting links
REGS <- GRNS2[,c("TF","chr","start","end","ID","Gene")]
save(REGS,file = paste0("SCENIC/",Cls,"_",Splt,"/GRNs/",Cls,"_",Splt,"_TF2R2G.RData"))

q()

## This is to merge the GRN from each split and find the GRN connections that appear in both
## In new run I have removed TF grom its GRNs else TF expression can drive signal
GRNa <- read.delim(paste0("SCENIC/",Cls,"_A/GRNs/",Cls,"_A_eRegulons_fix.txt"))
GRNb <- read.delim(paste0("SCENIC/",Cls,"_B/GRNs/",Cls,"_B_eRegulons_fix.txt"))

tfs <- intersect(unique(GRNa$TF),unique(GRNb$TF))
## the way to cbine the GRN is for each TF identify the targets that appear 
## in both. for this connection identify regions that appear in both.
## The TF to target importance and correltion are from same data  
GRNs <- c()
for(x in tfs){
  
  #extract connections for a TF
  delTF1 <- GRNa[GRNa$TF==x,]
  delTF2 <- GRNb[GRNb$TF==x,]
  
  #intersection of targets
  tar <- intersect(delTF1$Gene,delTF2$Gene)
  ## remove TF from its targets
  tar <- tar[!(tar %in%x)]
  GRN_sel <- c()
  for(y in tar){
    
    ## extracting connections, merging, selecting those that are duplicated
    ## chosing max values
    deltar1 <- delTF1[delTF1$Gene==y,]
    deltar2 <- delTF2[delTF2$Gene==y,]
    deltar <- rbind(deltar1,deltar2)
   
    deltar$X <- paste0(deltar$TF,"_",deltar$Region,"_",deltar$Gene)
    deltar <- deltar[duplicated(deltar$X),]
    deltar$X <- NULL
    GRN_sel <- rbind(GRN_sel,deltar)
    rm(deltar1,deltar2,deltar)
  }
  rm(y)
  GRNs <- rbind(GRNs,GRN_sel)
  rm(GRN_sel,delTF1,delTF2,tar)
  
  
}
write.table(GRNs,paste0("SCENIC/",Cls,"_eRegulons_fix.txt"),
            quote = FALSE, row.names = FALSE, sep ="\t")

##getting TF-GRN
def <- GRNs[,c("TF","Gene")]
tfs <- unique(def$TF)
regs <- list()
for(x in tfs){
  defx <-def[def$TF==x,]
  genes <- unique(defx$Gene)
  ## remove TF from targets
  genes <- genes[!(genes %in%x)]
  ## add TF in targerts
  #genes <-  unique(c(genes,x))
  if(length(genes)>19){
    regs[[x]] <-genes
  }
  
  rm(defx,genes)
}
rm(x)
sort(lengths(regs),decreasing = TRUE)

save(regs, file = paste0("SCENIC/",Cls,"_gmt.RData"))
q()

