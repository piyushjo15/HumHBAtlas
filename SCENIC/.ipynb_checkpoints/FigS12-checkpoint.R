suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggrepel)
    library(Matrix)
})

setwd("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/")


####### plottng overlapS for cofactor and coactivator ###########
load("SCENIC/AllGRNs_v2.RData")
all_tfs[all_tfs$Class=="Glioblast","Class"] <- "Glioblasts"
all_tfs$GRN <- paste0(all_tfs$TF,"_",all_tfs$Class)
row.names(all_tfs) <- all_tfs$GRN


## 3. Co-factors-------------
## checking if co-factors define celltype specific activity.
## focusing on ARNT2
## First I will check if the TFs that find to the same CREs and then TFs that regulate the same targets
DIR_GRN <- "~/MBsnANA/HBana/SCENICana/SCENICOut/Class/"
DIR_ATAC <- "~/MBsnANA/HBana/ATACana/Outs/"

tf <- "ARNT2"
all_tfsx <- all_tfs[all_tfs$TF==tf,]
#sel_cls <- all_tfsx$Class
#sel_grns <- all_tfsx$GRN
#sel_cls <- readLines(paste0("Tab5/",tf,"_ordofclass.txt"))
if(tf=="NFIC"){
  sel_cls <- c("GC","Neurons_A","Neuroblasts_C","CBINT",
               "MixedGlia","MES_VAS")
  sel_grns <- paste0(tf,"_",sel_cls)
}else if(tf=="ARNT2"){
  sel_cls <- c("Neurons_A","Neurons_B","Neuroblasts_B","Neuroblasts_C",
               "PreO","Progenitors","MixedGlia","Glioblasts")
  sel_grns <- paste0(tf,"_",sel_cls)
}
cofac <- list() ##co-factor, binding to same CRE
coact <- list() ##regulating same targets
for(x in sel_cls){
  df <- read.delim(paste0(DIR_GRN,x,"/SCENIC/",x,"_eRegulons_fix_v2.txt"))
  df1 <- df[df$TF==tf,c("TF","Region","Gene")]
 
  cres_sel <- unique(df1$Region)
  df2 <- df[df$Region %in% cres_sel,c("TF","Region","Gene")]
  tfa <- unique(df2$TF)
  tfa <- tfa[tfa!=tf]
  cofac[[x]] <- tfa
  rm(df2,cres_sel,tfa)
  
  genes_sel <- unique(df1$Gene)
  df3 <- df[df$Gene %in% genes_sel,c("TF","Region","Gene")]
  tfb <- unique(df3$TF)
  tfb <- tfb[tfb!=tf]
  coact[[x]] <- tfb
  rm(df3,genes_sel,tfb)

  rm(df,df1)
}
lengths(cofac)
lengths(coact)



overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (round(intersection/minlen,4))
}


##vectorized jaccard
OS_cofac<- outer(sel_cls, sel_cls, FUN = Vectorize(function(i, j) {
  overlapS(cofac[[i]], cofac[[j]])
}))
rownames(OS_cofac) <- colnames(OS_cofac) <- sel_cls

OS_coact <- outer(sel_cls, sel_cls, FUN = Vectorize(function(i, j) {
  overlapS(coact[[i]], coact[[j]])
}))
rownames(OS_coact) <- colnames(OS_coact) <- sel_cls



## heatmap
mb <- ((seq(1:101)-1)/100)

mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)
library(pheatmap)
hp1 <- pheatmap(OS_cofac,
                cluster_rows = FALSE,cluster_cols = FALSE,
                #show_rownames = FALSE,angle=90,
                show_colnames = FALSE,
                breaks = mb, color = mc2, fontsize = 8)
pdf(paste0("SuppFigs/",tf,"_OScofac.pdf"), width = 6.5, height = 6,pointsize=10)
print(hp1)
dev.off()
hp2 <- pheatmap(OS_coact,
                cluster_rows = FALSE,cluster_cols = FALSE,
                #show_rownames = FALSE,angle=90,
                show_colnames = FALSE,
                breaks = mb, color = mc2, fontsize = 8)
pdf(paste0("SuppFigs/",tf,"_OScoact.pdf"), width = 6.5, height = 6,pointsize=10)
print(hp2)
dev.off()

########

































