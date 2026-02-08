##analyzing the predicted class for set of selected clusters

setwd("~/MBsnANA/HBana/ISM/Outs/")

df <- read.delim("FetalBrain/pred_FB_deephb_lstmv3_57MP_53.txt",header = FALSE)
df <- read.delim("FetalBrain/pred_FB_deephb_lstmv2_57MP_73.txt",header = FALSE)
df <- read.delim("FetalBrain/pred_FB_deephb_cnn_57MP_99.txt",header = FALSE)

head(df)

beds <- read.delim("FetalBrain/MarkerCREsbed.txt",row.names = 1)
cres <- row.names(beds) 
row.names(df) <- cres

mp_cls <- readLines("~/MBsnANA/HBana/ISM/Scripts/TopicMPs.txt")
colnames(df) <- mp_cls

thres_mp <- read.delim("OutfrmGPU/MP_threshold_summary_deephb_lstmv3_57MP_53.txt")
thres_mp <- read.delim("OutfrmGPU/MP_threshold_summary_deephb_lstmv2_57MP_73.txt")
thres_mp <- read.delim("OutfrmGPU/MP_threshold_summary_deephb_cnn_57MP_99.txt")

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
save(pred,beds, file = "FetalBrain/pred_deep_lstmv3_53.RData")

save(pred,beds, file = "FetalBrain/pred_deep_lstmv2_73.RData")

save(pred,beds, file = "FetalBrain/pred_deep_cnn_99.RData")

#
### process--------------

#load("FetalBrain/pred_deep_lstmv3_53.RData")
load("FetalBrain/pred_deep_lstmv2_73.RData")
load("FetalBrain/pred_deep_cnn_99.RData")

## all MPs 
all_mp <- unique(unlist(pred))
all_mp <- all_mp[!is.na(all_mp)]
cls <- names(pred)
ct_ord <-readLines("FetalBrain/Celltype_order.txt")
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

sel1 <- readLines("FetalBrain/MP_order.txt")
sel1 <- paste0("r", mp_ord[sel1,"NewMP"])
sel <- readLines("FetalBrain/Celltype_order.txt")
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
png("FetalBrain/Figures/LSTMv3.png",units = "in",width = 6,height = 8,res = 300)
print(hp)
dev.off()


hp <- pheatmap::pheatmap(t(mat[,sel]),
                   clustering_method = "ward.D2",
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   
                   color = mc,
                   #breaks = mb,
                   #border_color = NA,
                   border_color = "grey",
                   fontsize=8)


ct_ord <- hp$tree_row
ct_ord <- ct_ord$labels[ct_ord$order]
mp_ord <- hp$tree_col
mp_ord <- mp_ord$labels[mp_ord$order]
write(ct_ord, "FetalBrain/Celltype_order.txt")
write(mp_ord, "FetalBrain/MP_order.txt")


# Apply the function to each list element and store the results
counts_list <- lapply(pred, count_characters)

# Convert the list of counts to a dataframe
df_counts <- bind_rows(counts_list, .id = "List_Name") #Use bind_rows

# Replace NA counts with 0.
df_counts[is.na(df_counts)] <- 0

# Print the resulting dataframe
print(df_counts)

# Example of how to do this using purrr and dplyr
df_counts_alt <- my_list %>%
  enframe(name = "List_Name", value = "character_list") %>%
  mutate(counts = map(character_list, ~{
    x <- .x[!is.na(.x)]
    table(factor(x, levels = unique(x)))
  })) %>%
  select(-character_list) %>% #remove the character_list column
  unnest_wider(counts, names_sep = ".") %>%
  mutate_if(is.numeric, replace_na, 0) #replaces NA with 0

print(df_counts_alt)

##