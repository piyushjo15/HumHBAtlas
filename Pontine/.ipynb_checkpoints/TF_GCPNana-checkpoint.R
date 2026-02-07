## Finding TF driving GC verus Pontine lineage
suppressPackageStartupMessages({
  library(Matrix)
  library(ggsci)
  library(ggrepel)
  library(tidyverse)
  library(scran)
})
DIR <- "~/MBsnANA/Diemana/scepro/"
setwd(DIR)

### --------
load("pdPNdwn.RData")
pdPN <- pd
load("pdGCdwn.RData")
pdGC <- pd
pd <- rbind(pdPN,pdGC)

load("lgcountsHB.RData")
load("~/MBsnANA/HBana/SCENICana/SCENICOut/GC/GC_HVG_pyscenic.RData")
tfs1 <- tfs_sel
load("~/MBsnANA/HBana/SCENICana/SCENICOut/PN/PN_HVG_pyscenic.RData")
tfs2 <- tfs_sel
tfs <- unique(c(tfs1,tfs2))

mdt <- mdt[tfs,row.names(pd)]
marks <- findMarkers(mdt, groups=pd$Annotationlv2, direction="up",
                     test.type="binom",pval.type="some")
cls <- names(marks)
sigs <-c()
for(x in cls){
  del <- marks[[x]]
  del <- row.names(del)[1:10]
  sigs <- rbind(sigs,del)
  rm(del)
}
rm(x)
row.names(sigs) <- cls

save(sigs, mdt, pd, file = "TF_GCPN.RData")
##----------
load("TF_GCPN.RData")
cl_lv <- c("Neurons_B:PN","Neurons_A:Pontine_A:a","Neurons_A:Pontine_C:a","Neurons_A:Pontine_B:a", 
            "Neuroblasts_A:A1_MES","Glioblast:Progenitor_MES",                     
            "GCUBC:GCP_UBCP","GCUBC:GCP_EM","GCUBC:GCP_PN",
            "GCUBC:GC_Diff_EM" ,"GCUBC:GC_Diff_PN"  , 
            "GC:GC_def","GC:GC_mature" )

cand <- readLines("del_PNGCTFGRNsel.txt")
cand <- unique(cand)
mdt <-mdt[,row.names(pd)]

cl_lv <- unique(pd$Annotationlv2)
cand <- unique(c(sigs))
plot.data <- pd
plot.data$RNN_clv2 <- pd$Annotationlv2
plotbub(cand)


##to obtain order
load("HBAnnlv2psb_280425.RData")
psb <- logcounts(del)
psb <- psb[cand,cl_lv]
hp <- pheatmap::pheatmap(psb,
                         scale = "row",
                         clustering_method = "ward.D2")
dela <- hp$tree_col 
cl_lv <- dela$labels[dela$order]
delb <- hp$tree_row 
cand <- delb$labels[delb$order]

######### plot functions----------
h <- length(cand)*0.3+2
w <- length(cl_lv)*0.3+2
#scaling 0 to 1
range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}


#new plot function
plotbub <- function(genes){
  exp_mat <- mdt[genes,]
  meta <- plot.data %>% 
    select(RNN_clv2)
  meta <- bind_cols(meta, as.data.frame(t(as.matrix(exp_mat))))
  head(meta) #View the first few lies
  meta <- pivot_longer(meta, -RNN_clv2, names_to="Gene", values_to="Expression")
  head(meta)
  metas <- meta %>%
    group_by(RNN_clv2, Gene) %>%
    summarise(Avg = mean(Expression),
              Pct = sum(Expression > 0) / length(Expression) * 100)
  metas2 <- metas %>%
    group_by(Gene) %>%
    mutate(ScaleAvg = range_01(Avg)) %>%
    ungroup()
  #head(metas)
  metas <- data.frame(metas2)
  metas$RNN_clv2 <- factor(metas$RNN_clv2, levels = cl_lv)
  metas$Gene <- factor(metas$Gene, levels = genes)
  p <- ggplot(metas, aes(x=RNN_clv2,y=Gene, size=Pct))+
    geom_point(aes(color=ScaleAvg))+
    scale_size("% detected", range = c(0,10)) +
    #scale_color_gradient(low = "white", high = "blue")+
    scale_color_viridis_c(direction = -1,option = "F" , 
                          name = "Scaled\nexpression")+
    theme_light()+
    theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
          axis.ticks = element_line(colour = 'black',linewidth=0.5),
          axis.text.x = element_text(face="bold",colour = "black",size=8,
                                     angle = 45,vjust=1, hjust=1),
          axis.text.y = element_text(face="bold",colour = "black",size=8),
          axis.title=element_blank())
  return(p)
  
}


########### TF-GRN analysis ############

load("AUC_MS/pdPNGC_PNGCGRNs_AUC.RData")
##downsample---
ncells <- 1000
cl_lv <- unique(pd$Annotationlv2)
cells <- c()
for(x in cl_lv){
  del_cells <- row.names(pd)[pd$Annotationlv2==x]
  if(length(del_cells)>1000){
    del_cells <- sample(del_cells,1000,replace = FALSE)
  }
  cells <- c(cells,del_cells)
  rm(del_cells)
}
pd <- pd[cells,]

## 2. New AUC bubble plot -----------
all_cells <- row.names(pd)
table(row.names(scr) %in%all_cells)
scr <- scr[all_cells,]

auc_rank_mat <- apply(scr, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC
table(all_cells==row.names(auc_rank_mat))
# Set high-AUC threshold (global or per-signature)
auc_co <- 0.95
mps <- colnames(scr)
# Create result list
res_list <- list()

for (x in mps) {
  auc_vec <- auc_rank_mat[, x]
  high_vec <- auc_vec > auc_co
  
  ## create a datafrom
  df <- data.frame(cell = rownames(auc_rank_mat), 
                   high_auc = high_vec, 
                   cluster = pd$Annotationlv2)
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
      TFGRN = x
    )
  
  res_list[[x]] <- res_sig
  rm(auc_vec,high_vec,all_high,res_sig)
}

# Combine and adjust p-values
merge_res <- bind_rows(res_list) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "BH"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-10))

merge_res <- merge_res %>%
  group_by(TFGRN) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 5), 0)) %>%
  ungroup()
ord_cl <- c("Neurons_B:PN","Neurons_A:Pontine_A:a","Neurons_A:Pontine_C:a","Neurons_A:Pontine_B:a", 
            "Neuroblasts_A:A1_MES","Glioblast:Progenitor_MES",                     
            "GCUBC:GCP_UBCP","GCUBC:GCP_EM","GCUBC:GCP_PN",
            "GCUBC:GC_Diff_EM" ,"GCUBC:GC_Diff_PN"  , 
            "GC:GC_def","GC:GC_mature" )
sel <- ord_sel <- readLines("PNGCTFGRNsel.txt")
merge_res <- data.frame(merge_res)
merge_res_sel <- merge_res[merge_res$TFGRN %in% sel,]
merge_res_sel$cluster <- factor(merge_res_sel$cluster,levels = ord_cl)
merge_res_sel$TFGRN <- factor(merge_res_sel$TFGRN,levels = ord_sel)

p <- ggplot(merge_res_sel, aes(x = cluster, y = TFGRN)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 95% AUC", range = c(0, 8)) +
  scale_color_gradient2(low = "white",mid = "grey87", high = "black",midpoint = 2.5,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        axis.text.x = element_text(face="bold",colour = "black",size=8,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title=element_blank())+
  labs(x = "MetaCREs", y = "Cluster",title = "Enrichment of TFGRN program per Cluster")
pdf("Figures/TFGRN_PNGC_prop.pdf", width = 6, height = 7, pointsize = 10)
print(p)
dev.off()





## marker analysis----------
pd$Group <- "Pontine"
cl_lv <- unique(pd$Annotationlv2)
sel <- grep("GC",cl_lv,value = TRUE)
pd[pd$Annotationlv2%in%sel,"Group"] <- "GC"

table(pd$Annotationlv2,pd$Group)
pd$Group2 <- "Imm"
sel <- c("GC","Neurons_A","Neurons_B")
pd[pd$Annotationlv1%in%sel,"Group"] <- "Mature"
pd$Group3 <- paste0(pd$Group,"_",pd$Group2)

scrx <- scr
pdx <- pd
## subset
pdx <- pd[pd$Group=="Pontine",]
scrx <- scr[row.names(pdx),]

marks <- scran::findMarkers(t(scrx), groups=pdx$Annotationlv2, direction="up",
                     test.type="wilcox",pval.type="all")
cls <- names(marks)
sigs <-c()
for(x in cls){
  del <- marks[[x]]
  del <- row.names(del)[1:5]
  sigs <- rbind(sigs,del)
  rm(del)
}
rm(x)
row.names(sigs) <- cls
sel <- unique(c(sigs))

### module score enrichment plots---------

zscore_c <- function(x){
  sdx <- sd(x)
  mx <- mean(x)
  ot <- boxplot.stats(x)$out
  q <- boxplot.stats(x)$stats
  keep <- (x < q[1])
  x[keep] <- q[1]
  keep <- (x > q[5])
  x[keep] <- q[5]
  x2 <- x-mx
  z <- x2/sdx
  return(z)
}
ms2x <- apply(scr, 2, zscore_c)
msd2 <- c()

cl_lv <- unique(pd$Annotationlv2)
for(x in cl_lv){
  keep1 <- pd$Annotationlv2==x
  msx <- ms2x[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]
msd3 <- apply(msd2,2,zscore_c)
msd3[1:4,1:4]

mc <- colorRampPalette(c("deepskyblue4","grey93","magenta3"))(101)
mb <- c(seq(-2, 0, length.out=ceiling(100/2) + 1),
        seq(2/100, 2, length.out=floor(100/2)))


hp <- pheatmap::pheatmap(msd3[,sel], 
                         #clustering_distance_cols = "correlation",
                         #clustering_distance_rows = "correlation",
                         #cluster_cols = FALSE,
                         color=mc, breaks=mb,
                         border_color = NA,
                         fontsize = 8)

dela <- hp$tree_col

ord_sel <- dela$labels[dela$order]
write(ord_sel,"del.txt")

pdf("Figures/MarkerTFGRNGCPN.pdf", width = 6, height = 6, pointsize = 10)
print(hp)
dev.off()

### distinct TF-GRNs in PN and GC lineage------
seta <- GSEA[["PN:ZNF483"]]
setb <- GSEA[["GC:ZNF483"]]

com <- intersect(seta,setb)
df <- data.frame(Genes=c(seta,setb),Group=c(rep("PN",length(seta)),rep("GC",length(setb))))
df <- df[!df$Genes%in%com,]
genes <- row.names(df) <- df$Genes

mdtx <- mdt[genes,]
ms2x <-t(apply(mdtx,1,zscore_c))


## extract data and pseudobulks
load("TF_GCPN.RData")
mdt <- del
mdtx <- mdt[,row.names(pd)]
rm(del)


cl_lv <- c("Neurons_B:PN","Neurons_A:Pontine_A:a",
            "GC:GC_def","GC:GC_mature" )

pdx <- pd[pd$Annotationlv2 %in% cl_lv,]
msd <- c()

for(x in cl_lv){
  cells <- row.names(pdx)[pdx$Annotationlv2==x]
  del <- rowMeans(ms2x[,cells])
  msd <- cbind(msd,del)
}
colnames(msd) <- cl_lv
msd[1:2,]
head(df)
df <- cbind(df,msd)
msd2 <- reshape2::melt(df)
head(msd2)
msd2$Var2 <- factor(msd2$variable,levels = cl_lv)
dim(msd2)

ggplot(msd2,aes(x=variable,y=value,fill = Group))+
  geom_boxplot()+
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        axis.text.x = element_text(face="bold",colour = "black",size=10),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title=element_blank())

dim(pdx)

cls <- c("GC","Neurons_A","Neurons_B")
pdy <- c()
for(y in cls){
  del <- pdx[pdx$Annotationlv1==y,]
  rn <- row.names(del)
  del2 <- del[sample(rn,1500),]
  pdy <- rbind(pdy,del2)
  rm(del,rn,del2)
}
rm(y)
pdx <- pdy
cells <- row.names(pdx)

## split into 50 cell blocks
arr <- sample(seq(1:length(cells)))
sizes <- rep(floor(length(cells)/ 50), 50)
remainder <- length(cells) %% 50

indices <- split(arr, rep(1:50, times = sizes))
pdx$Split <- "ND"
for(i in 1:length(indices)){
  pdx[indices[[i]],"Split"] <- paste0("a",i)
}
pdx$PSB <- paste0(pdx$Annotationlv1,"_",pdx$Split)

ms2x <-t(apply(mdt[genes,],1,zscore_c))
ms2x[1:4,1:4]
msd2 <- c()
cl_lv <- unique(pdx$PSB)
for(x in cl_lv){
  keep1 <- pdx$PSB==x
  msx <- ms2x[,keep1]
  del <- rowMeans(msx,na.rm = TRUE)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]


gene_grp <- factor(df$Group,levels = c("GC","PN")) ## row names of mp_cre_sel should match row-name of mat


ann <- pdx[!duplicated(pdx$PSB),c("Annotationlv1","PSB")] 
row.names(ann) <- ann$PSB
ann <- ann[row.names(msd2),] ### this is very important to check the colname order matches

cl_grp <- factor(ann$Annotationlv1,levels = cls)

col_class <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
value_class <- col_class[cls,"Color"]
names(value_class) <- cls

col_ha <- columnAnnotation(Group = cl_grp)
Heatmap(t(msd2),
              name = "Motif Enirchment",
              clustering_method_rows = "average",
              show_column_names = FALSE,
              show_row_names = FALSE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_split = gene_grp,
              column_split = cl_grp,top_annotatio =col_ha,
              cluster_row_slices = FALSE
        )
pdf("Figures/selMotifEnrichPreOOligo.pdf", width = 6, height = 4, pointsize = 1, useDingbats = FALSE)
draw(ht)
dev.off()
#### differentially expressed genes among the TF-GRN members------------

DIR_GRN <- "~/MBsnANA/HBana/SCENICana/SCENICOut/"
# TFGRNs
load(paste0(DIR_GRN,"PN/metadata/PN_gmt_v2.RData"))
PN_regs <- regs
names(PN_regs) <- paste0("PN:",names(PN_regs))
load(paste0(DIR_GRN,"GC/metadata/GC_gmt_v2.RData"))
GC_regs <- regs
names(GC_regs) <- paste0("GC:",names(GC_regs))
regs <- do.call(c,list(PN_regs,GC_regs))
GSEA <- regs
genes <- unique(unlist(GSEA))
rm(GSEA,regs)

load("TF_GCPN.RData")
load("lgcountsHB.RData")

pd$Group <- "Pontine"
cl_lv <- unique(pd$Annotationlv2)
sel <- grep("GC",cl_lv,value = TRUE)
pd[pd$Annotationlv2%in%sel,"Group"] <- "GC"
table(pd$Annotationlv2,pd$Group)
cl_lv <- c("Neuroblasts_A:A1_MES","Glioblast:Progenitor_MES",                     
           "GCUBC:GCP_UBCP","GCUBC:GCP_EM","GCUBC:GCP_PN" )
pdx <- pd[pd$Annotationlv2%in%cl_lv,]

marks <- scran::findMarkers(mdt[genes,row.names(pdx)], groups=pdx$Group, 
                            test.type="binom",pval.type="all")
df <- data.frame(marks$GC)
head(df)
tail(df)
seta <- row.names(df)[df$logFC.Pontine>1 & df$FDR<0.001]
setb <- row.names(df)[df$logFC.Pontine<1 & df$FDR<0.001]
write(seta,"del.txt")

pid_gene_sets = msigdbr(species = "Homo sapiens", collection = "H")

head(pid_gene_sets)
msigdbr_t2g = pid_gene_sets %>% dplyr::distinct(gs_exact_source,gene_symbol) %>% as.data.frame()
ks <- enricher(gene = seta, TERM2GENE = msigdbr_t2g)
res <- ks@result
res <- res[,c("ID","GeneRatio","p.adjust","geneID")]
res <- res[res$p.adjust<0.1,]






###