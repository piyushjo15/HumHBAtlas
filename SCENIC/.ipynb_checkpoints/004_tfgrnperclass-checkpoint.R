## This script process TF-GRNs per class to obtain enrichment values per TF
suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(AUCell)
  library(scran)
  library(pheatmap)
})

cls <- commandArgs(trailingOnly = TRUE)


setwd("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/")

range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}

## 1. load data ---------
load("SCENIC/AllGRNs_v2.RData")
all_tfs[all_tfs$Class=="Glioblast","Class"] <- "Glioblasts"
all_tfs$GRN <- paste0(all_tfs$TF,"_",all_tfs$Class)
row.names(all_tfs) <- all_tfs$GRN

all_tfs <- all_tfs[all_tfs$Class==cls,]
tfs <- all_tfs$TF
sel_grns <- all_tfs$GRN
write(tfs, file = paste0("Tab4/",cls,"_tfs.txt"))


load("SCENIC/pdHB_AllGRNs_MSz.RData")
kk <- colnames(scr_com)
colnames(scr_com) <- gsub("_Glioblast","_Glioblasts",kk,fixed = TRUE)
colnames(scr_com_z) <- gsub("_Glioblast","_Glioblasts",kk,fixed = TRUE)

load("plotdataRNA4Nilsv2.RData")
table(row.names(plot.data)==row.names(scr_com))
plot.data <- plot.data[plot.data$Class==cls,]
scr_com <- scr_com[row.names(plot.data),sel_grns]
colnames(scr_com) <- tfs
scr_com_z <- scr_com_z[row.names(plot.data),sel_grns]
colnames(scr_com_z) <- tfs


## 2. find top marker tf-grns per cluster---------
marks <- findMarkers(t(scr_com),groups=plot.data$Cluster, direction="up",
                     test.type="wilcox",pval.type="some")
cl <- names(marks)
sel1 <- c()
for(x in cl){
  del <- marks[[x]]
  del <- row.names(del)[1:2]
  sel1 <- rbind(sel1,del)
  rm(del)
}
rm(x)
sel1 <- unique(sel1)

marks <- findMarkers(t(scr_com),groups=plot.data$Cluster, direction="up",
                     test.type="wilcox",pval.type="all")
cl <- names(marks)
sel2 <- c()
for(x in cl){
  del <- marks[[x]]
  del <- row.names(del)[1:2]
  sel2 <- rbind(sel2,del)
  rm(del)
}
rm(x)
sel2 <- unique(sel2)

sel <- unique(c(sel1,sel2))

rm(marks,sel1,sel2)

## find order using enrichment
cl_lv <- unique(plot.data$Cluster)

msd2 <- c()

for(x in cl_lv){
  keep1 <- plot.data$Cluster==x
  msx <- scr_com_z[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
#msd2[1:4,1:4]


hp <- pheatmap(msd2[,sel], 
               clustering_method = "ward.D2",
               #clustering_distance_cols = "correlation",
               #clustering_distance_rows = "correlation",
               #cluster_rows = FALSE,
               silent = TRUE)
ord <- hp$tree_col
ord <- ord$labels[ord$order]
write(ord, file = paste0("Tab4/",cls,"_top2tfgrnpercluster.txt"))

rm(marks,all_tfs,all_grns)
## 2. Per cluster -----------

##max cells top 1000 cells per cluster
tb <- data.frame(table(plot.data$Cluster))
cl_lv <-row.names(tb)<- as.character(tb$Var1)

sel_cells <- c()
for(x in cl_lv){
  del_cell <- row.names(plot.data)[plot.data$Cluster==x]
  ll <- tb[x,2]
  if(ll>1000){
    del_cell <- sample(del_cell,1000)
  }
  sel_cells <- c(sel_cells,del_cell)
  rm(del_cell,ll)
}

pd <- plot.data[sel_cells,]
scr <- scr_com[sel_cells,]

## to a rank matrix first
auc_rank_mat <- apply(scr, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC

# Set high-AUC threshold (global or per-signature)
auc_co <- 0.90

# Create result list
res_list <- list()

for (x in tfs) {
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
      TFGRN = x
    )
  
  res_list[[x]] <- res_sig
  rm(auc_vec,high_vec,all_high,res_sig,df)
}

# Combine and adjust p-values
merge_res <- bind_rows(res_list) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p, method = "bonferroni"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-5))

merge_res <- merge_res %>%
  group_by(TFGRN) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 4), 0)) %>%
  ungroup()

merge_res <- data.frame(merge_res)

write.table(merge_res, file = paste0("Tab4/",cls,"_tfgrnenrpercluster.txt"),sep = "\t",quote = FALSE,row.names = FALSE)

rm(tb,merge_res,auc_rank_mat,res_list,auc_co,scr)
rm(sel_cells,pd)

## 3. Per Subcluster -----------

##max cells top 200 cells per subcluster
tb <- data.frame(table(plot.data$Subcluster))
cl_lv <-row.names(tb)<- as.character(tb$Var1)

sel_cells <- c()
for(x in cl_lv){
  del_cell <- row.names(plot.data)[plot.data$Subcluster==x]
  ll <- tb[x,2]
  if(ll>200){
    del_cell <- sample(del_cell,200)
  }
  sel_cells <- c(sel_cells,del_cell)
  rm(del_cell,ll)
}

pd <- plot.data[sel_cells,]
scr <- scr_com[sel_cells,]

## to a rank matrix first
auc_rank_mat <- apply(scr, 2, function(x){ return(rank(x,ties.method="first"))})
auc_rank_mat <- auc_rank_mat/dim(auc_rank_mat)[1] ## this converts rank to AUC

# Set high-AUC threshold (global or per-signature)
auc_co <- 0.90

# Create result list
res_list <- list()

for (x in tfs) {
  auc_vec <- auc_rank_mat[, x]
  high_vec <- auc_vec > auc_co
  
  ## create a datafrom
  df <- data.frame(cell = rownames(auc_rank_mat), 
                   high_auc = high_vec, 
                   cluster = pd$Subcluster)
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
  mutate(fdr = p.adjust(p, method = "bonferroni"))%>%
  mutate(log2OR = log2(odds_ratio + 1e-5))

merge_res <- merge_res %>%
  group_by(TFGRN) %>%
  mutate(log2OR_c =pmax(pmin(log2OR, 4), 0)) %>%
  ungroup()

merge_res <- data.frame(merge_res)

write.table(merge_res, file = paste0("Tab4/",cls,"_tfgrnenrpersubcluster.txt"),sep = "\t",quote = FALSE,row.names = FALSE)

q()
## 4. plot example------------
cls <- "PreO"

## load cluster and subcluster order
load("Tab1/OrderofClusters.RData")
ord_cl <- ord_cluster[[cls]]
ord_subcl <- ord_subcluster[[cls]]

## per cluster
merge_res <- read.delim(paste0("Tab4/",cls,"_tfgrnenrpercluster.txt"),header = TRUE)
sel <- readLines(paste0("Tab4/",cls,"_top2tfgrnpercluster.txt"))



merge_res_sel <- merge_res[merge_res$TFGRN %in% sel,]
merge_res_sel$cluster <- factor(merge_res_sel$cluster,levels = ord_cl)
merge_res_sel$TFGRN <- factor(merge_res_sel$TFGRN,levels = sel)

p <- ggplot(merge_res_sel, aes(x = cluster, y = TFGRN)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 8)) +
  scale_color_gradient2(low = "white",mid = "grey87", high = "black",midpoint = 2,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        axis.text.x = element_text(face="bold",colour = "black",size=10,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank())+
  labs(
    x = "TF-GRN",
    y = "Cluster",
    title = "Enrichment of TFGRN program per Cluster",
  )

rm(merge_res_sel,merge_res,p)
## per subcluster
merge_res2 <- read.delim(paste0("Tab4/",cls,"_tfgrnenrpersubcluster.txt"),header = TRUE)

merge_res_sel <- merge_res2[merge_res2$TFGRN %in% sel,]
merge_res_sel$cluster <- factor(merge_res_sel$cluster,levels = ord_subcl)
merge_res_sel$TFGRN <- factor(merge_res_sel$TFGRN,levels = sel)

p <- ggplot(merge_res_sel, aes(x = cluster, y = TFGRN)) +
  geom_point(aes(size = prop, color = log2OR_c)) +
  scale_size_continuous(name = "Proportion\n in top 90% AUC", range = c(0, 8)) +
  scale_color_gradient2(low = "white",mid = "grey87", high = "black",midpoint = 2,
                        name = "Log2 OR")+
  #theme_minimal(base_size = 12) +
  theme_light()+
  theme(axis.line = element_line(colour = 'black',linewidth = 0.5),
        axis.ticks = element_line(colour = 'black',linewidth=0.5),
        axis.text.x = element_text(face="bold",colour = "black",size=10,
                                   angle = 45,hjust=1),
        axis.text.y = element_text(face="bold",colour = "black",size=10),
        axis.title=element_blank())+
  labs(
    x = "TF-GRN",
    y = "Cluster",
    title = "Enrichment of TFGRN program per Cluster",
  )

