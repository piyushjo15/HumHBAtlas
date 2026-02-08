
###fig S11 C
load("~/MBsnANA/HBana/ISM/Outs/FetalBrain/pred_deep_lstmv2_73.RData")

## all MPs 
all_mp <- unique(unlist(pred))
all_mp <- all_mp[!is.na(all_mp)]
cls <- names(pred)
ct_ord <-readLines("~/MBsnANA/HBana/ISM/Outs/FetalBrain/Celltype_order.txt")
df_counts <- c()

for(x in cls){
  del_pred <- pred[[x]]
  del_pred <- del_pred[!is.na(del_pred)]
  counts <- table(factor(del_pred, levels = all_mp))
  counts <- counts/20
  df_counts <- cbind(df_counts,counts)
  rm(counts)
}
colnames(df_counts) <- cls
df_counts[1:4,1:4]

rs <- apply(df_counts, 1, max)
keep <- rs>10
df_counts_sel <- df_counts[keep,]
mat <- as.matrix(df_counts_sel)
mc <- colorRampPalette(c("white","grey50","black"))(101)
#mb <- seq(1:101)-1

mat[1:4,1:4]
mat[1:4,1:4]
mp_ord <- read.delim("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/OrderofMPsforplots.txt",row.names = 1)
cn <- rownames(mat)
rownames(mat) <- mp_ord[cn,"NewMP"]
rownames(mat)<- paste0("r", rownames(mat))
mat[1:4,1:4]
## select columns
# cn <- colnames(mat)
# df <- read.delim("FetalBrain/Cluster_metadata.txt")
# df[1:4,1:4]
# df$X <- NULL
# df$Clusters <- paste0("Cl",df$Clusters)
# sel_class <- c("Oligo","Immune")
# sel1 <- df$ClusterName[df$Class %in% sel_class]
# sel2 <- grep("CB_GNP",df$ClusterName,value = TRUE)
# sel3 <- grep("Purk",df$ClusterName,value = TRUE)
# sel <- unique(c(sel1,sel2,sel3))
##select rMP

sel1 <- readLines("~/MBsnANA/HBana/ISM/Outs/FetalBrain/MP_order.txt")
sel1 <- paste0("r", mp_ord[sel1,"NewMP"])
sel <- readLines("~/MBsnANA/HBana/ISM/Outs/FetalBrain/Celltype_order.txt")
#select cluster
mc <- colorRampPalette(c("white","grey50","black"))(101)

hp <- pheatmap::pheatmap(t(mat[sel1,sel]),
                         #clustering_method = "ward.D2",
                         #clustering_distance_cols = "correlation",
                         #cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         color = mc,
                         border_color = "grey",
                         fontsize=8)


pdf("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/SuppFigs/FetBrainDeepHBpred.pdf",width = 6,height = 8,pointsize = 10)
print(hp)
dev.off()
##CORRELATION
load("~/MBsnANA/HBana/ATACana/Outs/CREs/MergedpdAUC_FB_CRE.RData")
mdt_fb <- mdt
load("~/MBsnANA/HBana/ATACana/Outs/CREs/MergedpdAUC_MPCRE.RData")
mdt_mp <- mdt2
dim(mdt_fb)
dim(mdt_mp)
rm(mdt,mdt2)
rn <- intersect(row.names(mdt_fb),row.names(mdt_mp))
mdt_fb <- mdt_fb[rn,]
mdt_mp <- mdt_mp[rn,]
mp_ord <- read.delim("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/OrderofMPsforplots.txt",row.names = 1)
cn <- colnames(mdt_mp)
colnames(mdt_mp) <- paste0("r", mp_ord[cn,"NewMP"])

sel1 <- paste0("rMP", c(33,53,54,52,55,24,27,32,20,37,39))
sel <- readLines("~/MBsnANA/HBana/ISM/Outs/FetalBrain/Celltype_order.txt")

mdt_mp <- mdt_mp[,sel1]
mdt_fb <- mdt_fb[,sel]

corx <- cor(mdt_fb,mdt_mp)
mc <- colorRampPalette(c("blue4","blue2","white","orangered2","orangered3"))(201)
mb <- ((seq(1:201)-1)/100)-1
hp <- pheatmap::pheatmap(corx[sel,sel1], 
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         clustering_method = "ward.D2",
                         color=mc, breaks=mb,
                         legend = FALSE,
                         border_color = "grey",
                         fontsize = 12)

#### Fig S11B #############
setwd("~/MBsnANA/HBana/ISM/Outs/OutfrmGPU/")

## on all the topic CREs
cnn <- read.delim("ROC_PR_deephb_cnn_99.txt")
cnn$Model <- "CNN"
lstmv2 <- read.delim("ROC_PR_deephb_lstmv2_73.txt")
lstmv2$Model <- "lstmv2"

lstmv3 <- read.delim("ROC_PR_deephb_lstmv3_53.txt")
lstmv3$Model <- "lstmv3"

cnn_a <- read.delim("MP_threshold_summary_deephb_cnn_57MP_99.txt")
cnn <- cbind(cnn_a,cnn)
lstmv2_a <- read.delim("MP_threshold_summary_deephb_lstmv2_57MP_73.txt")
lstmv2 <- cbind(lstmv2_a,lstmv2)

lstmv3_a <- read.delim("MP_threshold_summary_deephb_lstmv3_57MP_53.txt")
lstmv3 <- cbind(lstmv3_a,lstmv3)


## only on test split-----
cnn <- read.delim("HBCRE_test/ROC_PR_deephb_cnn_99.txt")
cnn$Model <- "CNN"

lstmv2 <- read.delim("HBCRE_test/ROC_PR_deephb_lstmv2_73.txt")
lstmv2$Model <- "lstmv2"

lstmv3 <- read.delim("HBCRE_test/ROC_PR_deephb_lstmv3_53.txt")
lstmv3$Model <- "lstmv3"

cnn_a <- read.delim("HBCRE_test/MP_threshold_summary_deephb_cnn_57MP_99.txt")
cnn <- cbind(cnn_a,cnn)
lstmv2_a <- read.delim("HBCRE_test/MP_threshold_summary_deephb_lstmv2_57MP_73.txt")
lstmv2 <- cbind(lstmv2_a,lstmv2)

lstmv3_a <- read.delim("HBCRE_test/MP_threshold_summary_deephb_lstmv3_57MP_53.txt")
lstmv3 <- cbind(lstmv3_a,lstmv3)

mp_ord <- read.delim("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/OrderofMPsforplots.txt",row.names = 1)
head(mp_ord)
ord <- mp_ord$NewMP<- paste0("r", mp_ord$NewMP)


merge_df <- rbind(cnn,lstmv2,lstmv3)
merge_df$OldMP <- merge_df$MP
mps <- unique(merge_df$MP)
for(x in mps){
  merge_df[merge_df$OldMP==x,"MP"] <- mp_ord[x,"NewMP"]
}
merge_df$OldMP <- NULL
merge_df$X <- merge_df$Best.Threshold <- merge_df$Precision <- merge_df$Recall <- NULL
merge_df2 <- reshape2::melt(merge_df)
head(merge_df2)
tail(merge_df2)
merge_df2$Model <- factor(merge_df2$Model, levels = c("CNN","lstmv3","lstmv2"))
merge_df2$variable <- factor(merge_df2$variable, levels = c("PR","ROC","F1.Score"))
p <- ggplot(merge_df2,aes(x=variable,y=value,fill=Model))+
  geom_boxplot(color="blue")+
  scale_fill_manual(values=c("white","grey75","black"))+
  theme_bw()+
    theme(
          axis.ticks = element_line(colour = 'black',linewidth=.5),
          axis.text.x = element_text(face="bold",colour = "black",size=20),
          axis.text.y = element_text(face="bold",colour = "black",size=15),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          axis.title=element_blank())




## bar plot per metric per model
merge_df2$MP <- factor(merge_df2$MP,levels=rev(ord))

p <- ggplot(merge_df2,aes(y=MP,x=value,fill=Model))+
  geom_bar(stat = "identity",width = 0.6,color="blue",position = position_dodge())+
  scale_fill_manual(values=c("white","grey75","black"))+
  theme_bw()+
  theme(
    axis.ticks = element_line(colour = 'black',linewidth=.5),
    axis.text.x = element_text(colour = "black",size=14),
    axis.text.y = element_text(colour = "black",size=14),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.title=element_blank())+
  scale_x_continuous(breaks = c(0,0.5,1))+
  facet_wrap(~variable)


pdf("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/SuppFigs/FigS11b.pdf", width = 6, height = 16,pointsize = 10)
print(p)
dev.off()