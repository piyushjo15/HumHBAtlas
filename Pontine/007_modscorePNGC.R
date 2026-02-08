##This script calculates enrichment scores (using Seurat algorithm) for a gene-set
## for hindbrain data s

suppressPackageStartupMessages({
  library(Matrix)
  library(scater)
  library(scran)
  library(dplyr)
  library(tidyr)
  library(AUCell)
  library(BiocParallel)
})
addmodscore <- function()
## function defined previously, not shown here for space
##Directories 
DIR_RNA <- "~/RNA/"
DIR_GRN <- "~/SCENIC/"
DIR <- "~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/PNana/"
setwd(DIR)
# 
##1. load required data -------
### 1.1 load Gene-sets--------
# ## genesets
# load(paste0(DIR_RNA,"Reactome_Signaling_Pathways.RData"))
# ### remove ribosomal and mitochndiral genes
# load(paste0(DIR_RNA,"genesnoMTRibo.RData"))
# GSEAx <- GSEA
# GSEA <- list()
# cls <- names(GSEAx)
# for(x in cls){
#   cand <- GSEAx[[x]]
#   cand <- cand[cand %in% genes]
#   if(length(cand)>19){
#     GSEA[[x]] <- cand
#   }
#   rm(cand)
# }
# 
# # TFGRNs
# load(paste0(DIR_GRN,"PN/metadata/PN_gmt.RData"))
# PN_regs <- regs
# names(PN_regs) <- paste0("PN:",names(PN_regs))
# 
# load(paste0(DIR_GRN,"Class/GCUBC/SCENIC/GCUBC_gmt.RData"))
# GCUBC_regs <- regs
# names(GCUBC_regs) <- paste0("GCUBC:",names(GCUBC_regs))
# 
# regs <- do.call(c,list(PN_regs,GCUBC_regs))
# GSEA <- regs


## top 1000 DEG between pontine nuclei or granule cell lienage
load("DEG_PNGCtop100.RData")

## 1.2 load expression data -----------
## this is compiled log transformed expression data for the entire hindbrain
load("~/RNA/lgcountsHB.RData")

# ## obtain and subset metadata ----------
# # ##load metadata, combined 
load("~/RNA/plotdataHBv2.RData")

cl_lv <- c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
pd <- plot.data[plot.data$Cluster %in% cl_lv,]
cells <- row.names(pd)

mdt <- mdt[,cells]

## 2. module score----
scr <- addmodscore(mdt,features = GSEA)
colnames(scr) <- names(GSEA)
row.names(scr) <- row.names(pd)


#save(scr, pd, file = paste0("pdPNGCUBC_GRNs_MS.RData")) ## PN and GCUBC GRNss
#save(scr, pd, file = paste0("pdPNGC_RetSig_MS.RData")) ## Reactome Signalome
save(scr, pd, file = paste0("pdPNGC_DEGt100PNGC_MS.RData")) ## 

q()
##### Analysis ################
## 4. AUC score amrker analysis --------
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  }

load("pdPNGCUBC_GRNs_MS.RData")
#load("pdPNGC_RetSig_MS.RData")
#load(paste0(DIR_RNA,"Reactome_Signaling_Pathways.RData"))
rm(GSEA)
plot.data <- pd

##susbet for markers
cl_lv <- c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
pd <- plot.data[plot.data$Cluster %in% cl_lv,]
cells <- row.names(pd)
pd$Group <- "GC"
pd[pd$Cluster %in% c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES"),"Group"] <- "PN"
scr <- scr[row.names(pd),]

## marker by group
marks <- scran::findMarkers(t(scr),groups=pd$Group, direction="up",
                     test.type="wilcox",pval.type="some")
marks <- scran::findMarkers(t(scr),groups=plot.data$Cluster, direction="up",
                            test.type="wilcox",pval.type="some")
cl <- names(marks)
sig <- list()

for(x in cl){
  del <- marks[[x]]
  del <- del[del$FDR<1e-5 & del$summary.AUC>0.85,]
  del <- row.names(del)
  if(x %in% cl[1:2]){
    keep <- grep("PN:",del)
    if(length(keep)!=0){
      del <- del[-c(keep)]
    }
  }else{
    keep <- grep("GCUBC:",del)
    if(length(keep)!=0){
      del <- del[-c(keep)]
    }
  }
  sig[[x]] <-del
  rm(del)
}
rm(x)
lengths(sig)
sel <- unique(unlist(sig))

write(sel, file = "TF_GRNS_PN_GCUBC.txt")
## already hav downsamples to 1500 cells
## I would plot two GRNs, one full set of differential TFGRNS for Supp fig
## and one selected for main panel focusing on WNT and HOX

load("pdPNGCUBC_GRNs_MS.RData")

##downsample for enrichment plot
cells <- c()
for(x in cl_lv){
  keep <- row.names(pd)[pd$Cluster==x]
  if(length(keep)>1500){
    keep <- sample(keep,1500)
  }
  cells <- c(cells,keep)
  rm(keep)
}

pd <- pd[row.names(pd) %in% cells,]
table(pd$Cluster)
cells <- row.names(pd)
write(cells, file = "PN_GCUBC_cells.txt")
    
cl_lv <- c("Neuroblasts_I:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
pd <- pd[pd$Cluster %in% cl_lv,]
scr <- scr[row.names(pd),]
##for reactome and signalome
scr2 <- scr[,grn]

auc_rank_mat <- apply(scr2, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC
table(row.names(pd)==row.names(auc_rank_mat))
# Set high-AUC threshold (global or per-signature)
auc_co <- 0.90

# Create result list
res_list <- list()

for (x in grn) {
  auc_vec <- auc_rank_mat[, x]
  high_vec <- auc_vec > auc_co
  
  ## create a datafrom
  df <- data.frame(cell = rownames(auc_rank_mat), 
                   high_auc = high_vec, 
                   cluster = pd$Cluster)
  row.names(df) <- df$cell
  # Total cell above rank cut-off across clusters
  all_high <- sum(df$high_auc)
  ncells <- nrow(df)
  
  res_sig <- df %>%
    group_by(cluster) %>%
    summarise(
      tot = n(), ### total number of cells in a cluster
      a = sum(high_auc), # In cluster, high AUC
      b = tot - a,# In cluster, not high AUC
      c = all_high-a,# Outside cluster, high AUC
      d = (ncells- tot)-c,
      prop = (a / tot)*100,# % of above cut-off cells
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      odds_ratio = (a * d) / (b * c + 1e-10),
      p = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value,
      MP = x
    )
  
  res_list[[x]] <- res_sig
  rm(auc_vec,high_vec,ncells,all_high,res_sig)
}

# Combine and adjust p-values
merge_res <- bind_rows(res_list) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "BH"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-10))

##for sel
merge_res <- merge_res %>%
  group_by(MP) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 2.5), 0)) %>%
  ungroup()
head(merge_res)
merge_res <- data.frame(merge_res)

##selected GRNS  ----
ord_mp <- c("PN:TCF7L1","PN:TCF7L2","PN:ELF1","PN:HMG20A","PN:SREBF2",
                   "PN:HOXA3","PN:HOXB2","PN:HOXB3","PN:PAX6","GCUBC:PAX6")




merge_res$cluster <- factor(merge_res$cluster,levels = cl_lv)
merge_res$MP <- factor(merge_res$MP,levels = ord_mp)

p <- ggplot(merge_res, aes(x=cluster, y= MP)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 10)) +
  scale_color_gradient2(low = "white",mid = "grey87", high = "black",midpoint = 1.25,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        
        axis.text.x = element_text(colour = "black",size=12,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(colour = "black",size=12),
        axis.title=element_blank())



########## Reactome ##############

### AUC bubble for Ret Sig -----------
load("pdPNGC_RetSig_MS.RData")
load(paste0(DIR_RNA,"Reactome_Signaling_Pathways.RData"))
rm(GSEA)
##susbet for markers
cl_lv <- c("Neuroblasts_I:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
pd$Group <- "GC"
pd[pd$Cluster %in% c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES"),"Group"] <- "PN"
scr <- scr[row.names(pd),]

## marker by group
marks <- scran::findMarkers(t(scr),groups=pd$Group, direction="up",
                     test.type="wilcox",pval.type="some")
cls <- names(marks)
sig <- list()
for(x in cls){
  del <- marks[[x]]
    del <- del[del$summary.AUC>0.7 & del$FDR<1e-4,]
  del <- row.names(del)
  sig[[x]] <- del
  rm(del)
}
rm(x)

sel <- unique(unlist(sig))

del2 <- part2[sel,]
write.table(del2, file = "orderRetSigN.txt",sep="\t",quote = FALSE,row.names = FALSE)
write(sel, file = "RetSig_PN_GCUBC.txt")

##fixed for order and removed redundancies
del <- read.delim("orderRetSig.txt",row.names = 2)
sel <- row.names(del)
load("pdPNGC_RetSig_MS.RData")
cells <- readLines("PN_GCUBC_cells.txt")
pd <- pd[cells,]
scr <- scr[row.names(pd),]
scr2 <- scr[,sel]
ord_mp <- mps <- colnames(scr2) <- del[sel,"displayName"]
scr2[1:4,1:4]

auc_rank_mat <- apply(scr2, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC
table(row.names(plot.data)==row.names(auc_rank_mat))
# Set high-AUC threshold (global or per-signature)
auc_co <- 0.90

# Create result list
res_list <- list()

for (x in mps) {
  auc_vec <- auc_rank_mat[, x]
  high_vec <- auc_vec > auc_co
  
  ## create a datafrom
  df <- data.frame(cell = rownames(auc_rank_mat), 
                   high_auc = high_vec, 
                   cluster = pd$Cluster)
  row.names(df) <- df$cell
  # Total cell above rank cut-off across clusters
  all_high <- sum(df$high_auc)
  ncells <- nrow(df)
  
  res_sig <- df %>%
    group_by(cluster) %>%
    summarise(
      tot = n(), ### total number of cells in a cluster
      a = sum(high_auc), # In cluster, high AUC
      b = tot - a,# In cluster, not high AUC
      c = all_high-a,# Outside cluster, high AUC
      d = (ncells- tot)-c,
      prop = (a / tot)*100,# % of above cut-off cells
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      odds_ratio = (a * d) / (b * c + 1e-10),
      p = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value,
      MP = x
    )
  
  res_list[[x]] <- res_sig
  rm(auc_vec,high_vec,ncells,all_high,res_sig)
}

# Combine and adjust p-values
merge_res <- bind_rows(res_list) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "BH"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-10))

merge_res <- merge_res %>%
  group_by(MP) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 3), 0)) %>%
  ungroup()

head(merge_res)
merge_res <- data.frame(merge_res)

merge_res$cluster <- factor(merge_res$cluster,levels = cl_lv)
merge_res$MP <- factor(merge_res$MP,levels = sel)




p <- ggplot(merge_res, aes(x=cluster, y= MP)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 10)) +
  scale_color_gradient2(low = "white",mid = "grey70", high = "black",midpoint = 1.5,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        
        axis.text.x = element_text(colour = "black",size=12,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(colour = "black",size=12),
        axis.title=element_blank())+

