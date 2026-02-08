##analyzing the predicted class for set of selected clusters

setwd("~/DeepHB/Outs/")

df <- read.delim("FetalBrain/pred_FB_deephb_lstmv2_57MP_73.txt",header = FALSE)

head(df)

beds <- read.delim("FetalBrain/MarkerCREsbed.txt",row.names = 1)
cres <- row.names(beds) 
row.names(df) <- cres

mp_cls <- readLines("~/DeepHB/TopicMPs.txt")
colnames(df) <- mp_cls

thres_mp <- read.delim("MP_threshold_summary_deephb_lstmv2_57MP_73.txt")

row.names(thres_mp) <- thres_mp$MP
thres_mp <- thres_mp$Best.Threshold
names(thres_mp) <- mp_cls

beds$ID <- NA
## create a matrix of 0 and 1
df2 <- matrix(0,dim(df)[1],dim(df)[2])
colnames(df2)<- colnames(df)
row.names(df2)<- row.names(df)
for (i in 1:ncol(df)) {
  df2[, i] <- ifelse(df[, i] >= thres_mp[i], 1,0)
}

rs <- rowSums(df2)
# ## find which cres didn't get assigned any id and remove them
# keep <- rs==0
# cres_rem <- cres[keep]
## which CREs have only one prediction
keep <- rs==1
cres_1c<- cres[keep]
sel_ind <- max.col(df2[cres_1c,])

beds[cres_1c,"ID"] <-mp_cls[sel_ind]

## now for CREs with more than one prediction,
## find the prediction that passes the cut-off and is max
keep <- rs>1
cres_mc<- cres[keep]

## make below cut-off as 0
df2x <- df2[cres_mc,]
dfx <- df[cres_mc,]
keep <-df2x==0
dfx[keep]<- 0
sel_ind <- max.col(dfx)
beds[cres_mc,"ID"] <-mp_cls[sel_ind]



## identifying which class got which MP--------
fn <- list.files("FetalBrain/markerpeaks/")
fn <- gsub("_peaks.bed","",fn,fixed = TRUE)

pred <- list()

for(x in fn){
  del <- read.delim(paste0("FetalBrain/markerpeaks/",x,"_peaks.bed"),header = FALSE)
  del <- del$V4
  bedx <- beds[del,]
  pred[[x]] <- bedx$ID
  rm(del,bedx)
  
}

rm(x)

save(pred,beds, file = "FetalBrain/pred_deep_lstmv2_57MP_73.RData")

load("FetalBrain/pred_deep_lstmv2_57MP_73.RData")

## all MPs 
all_mp <- unique(unlist(pred))
all_mp <- all_mp[!is.na(all_mp)]
cls <- names(pred)
## select cell types
# cn <- colnames(mat)
# df <- read.delim("FetalBrain/Cluster_metadata.txt")
# df[1:4,1:4]
# df$X <- NULL
# df$Clusters <- paste0("Cl",df$Clusters)
# sel_class <- c("Oligo","Immune")
# sel1 <- df$ClusterName[df$Class %in% sel_class]
# sel2 <- grep("CB_GNP",df$ClusterName,value = TRUE)
# sel3 <- grep("Purk",df$ClusterName,value = TRUE)
# ct_ord <- unique(c(sel1,sel2,sel3))

ct_ord <- c("Neur_Purk_4","Neur_Purk_2" ,"Neur_Purk_3","Neur_Purk_5","Neur_Purk_1","Neur_Purk_6",
        "Mic_active_1","PVM", "Mic_2","Mic_3","Mic_dividing","Mic_5","Mic_6","Mic_active_2","Mic_1",
        "Immune_uncertain","Mic_4","Schwl_2", "Schwl_3","PreOPC_1","Schwl_1",
        "T_Cell_immature","Neur_CB_GNP_IPC_2","Neur_CB_GNP_IPC_4","Neur_CB_GNP_IPC_1",
        "Neur_CB_GNP_IPC_3","Neur_Purk_NBL_1","Neur_Purk_NBL_2","Neur_Purk_prog",
        "COP","OPC_6","OPC_2","OPC_1","OPC_3","PreOPC_3","PreOPC_4","OPC_4","OPC_5",
        "PreOPC_2")

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

mat[1:4,1:4]

mp_ord <- c("rMP33","rMP53","rMP54","rMP52","rMP55",
            "rMP24","rMP27","rMP32","rMP20","rMP37","rMP39")

#select cluster
mc <- colorRampPalette(c("white","grey50","black"))(101)

hp <- pheatmap::pheatmap(t(mat[mp_ord,ct_ord]),
                            cluster_rows = FALSE,
                         color = mc,
                         border_color = "grey",
                         fontsize=8)


