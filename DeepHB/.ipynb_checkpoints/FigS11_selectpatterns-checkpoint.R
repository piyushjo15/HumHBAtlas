library(tidyverse)

setwd("~/DeepHB/Outs/")

## select patterns based on some filter
pat_names <- readLines("~/MBsnANA/HBana/ISM/Outs/OutfrmGPU/all_pat_names.txt")
seq_ids <- paste0("seqid_",sprintf("%03d", 1:165))
pat_names <- data.frame(X=pat_names)
pat_names <- pat_names %>% separate(X,c("A","B"))
pat_names <- pat_names$B
pat_list <- str_split(pat_names,"")
keep <- lengths(pat_list)>3
table(keep)
seq_ids <- seq_ids[keep]
##old to new mp conver
mpconv <- read.delim("MPidconversion.txt")
oldmp <- row.names(mpconv) <- mpconv$MPs
new_mp <- mpconv$NewMP

## plotting full motif heatmap

### indentify imporant motifs from seqcount
load( "~/MBsnANA/HBana/ATACana/Outs/extrafiles/Seqletlogcounts_allPat.RData")

pattern_matrix <- pattern_matrix[oldmp,seq_ids]
row.names(pattern_matrix) <- new_mp


## remvoe neg values
pattern_matrix[pattern_matrix<0] <- 0
## plotting selected motif heatmap
rs <- colSums(pattern_matrix>0)
keep <- rs>0
table(keep)
sel_seq_ids <- colnames(pattern_matrix)[keep]

pat_mat <- pattern_matrix[,sel_seq_ids]
dim(pat_mat)

pat_mat2 <- matrix(0,nrow(pat_mat),ncol(pat_mat))
pat_mat2[pat_mat>0]<-1
row.names(pat_mat2) <- row.names(pat_mat)
colnames(pat_mat2) <- colnames(pat_mat)

cs <- data.frame(Freq= colSums(pat_mat2))
cs$Seqid <- row.names(cs)
cs <- cs[order(cs$Freq,decreasing = TRUE),]
cs$Seqid <- factor(cs$Seqid, levels=cs$Seqid)
head(cs)
tb <- data.frame(table(cs$Freq))
head(tb)

##subjectively divide into classes
cs$Class <- "A"
cs[cs$Freq>2 & cs$Freq<9,"Class"] <- "B"
cs[cs$Freq>=9,"Class"] <- "C"
table(cs$Class)

#plot
tb$X <- as.numeric(tb$Var)
tb$Class <- "A"
tb[tb$X>2 & tb$X<9,"Class"] <- "B"
tb[tb$X>=9,"Class"] <- "C"
cx <- c("beige","bisque1","bisque4")
names(cx) <- c("A","B","C")
p <- ggplot(tb, aes(x=Var1,y=Freq, fill=Class))+
scale_fill_manual(values=cx)+
geom_bar(stat="identity",color="black")+
theme_bw()

pdf("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/SuppFigs/FigS11E.pdf", width = 8, height = 4, pointsize = 10)
print(p)
dev.off()
##hetmap
ann <- cs
ann$Seqid <- ann$Freq <- NULL
ann_col <- list(Class=cx)

mc <- colorRampPalette(c("white","grey70","black"))(100)
hp<-pheatmap(t(pat_mat),
             annotation_row =  ann,annotation_colors=ann_col,
             clustering_method = "centroid",
             color=mc,
             cluster_cols=FALSE,border_color = NA)

seqidord <- hp$tree_row
seqidord <- seqidord$labels[seqidord$order]
head(seqidord)


pdf("~/MBsnANA/HBMS_imp/HB_copy_of_imp_data/SuppFigs/FigS11D.pdf", width = 12, height = 26, pointsize = 10)
print(hp)
dev.off()

## also identifying top pattern per MP
## finding top 3 motifs based on importance per MP
## also identifying top pattern per MP
## finding top 3 motifs based on importance per MP
seq_id_list <- list()
seq_id_mat <- c()
seq_ids_sel <- c()
for(x in new_mp){
  del <- names(sort(pattern_matrix[x,],decreasing = TRUE))[1:2]
  seq_id_list[[x]] <- del
  seq_id_mat <- rbind(seq_id_mat,del)
  seq_ids_sel <- c(seq_ids_sel,del[1])
  rm(del)
}
rm(x)
row.names(seq_id_mat) <- new_mp
seq_ids_sel <- unique(seq_ids_sel) ## 65 in top 3, 43 in top 2, 28 in top 1
sel_seq_ids <- unique(unlist(seq_id_list))
length(sel_seq_ids)

row.names(seq_id_mat) <- new_mp
write.table(seq_id_mat, file = "top3seqidperMP.txt",sep="\t",quote = FALSE)
write.table(seq_id_mat, file = "top2seqidperMP.txt",sep="\t",quote = FALSE)

save(seq_id_list, file = "topseqidsperMP.RData")
write(seq_ids_sel, file = "topseqidperMP.txt")



write(sel_seq_ids, file ="~/MBsnANA/HBana/ISM/Outs/OutfrmGPU/selectedseqidsforplot.txt" )
pat_mat <- pattern_matrix[,sel_seq_ids]

range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
pat_mat3 <- apply(pat_mat,2,range_01)

hp<- pheatmap::pheatmap(t(pat_mat3[,sel_seq_ids]),
                        ,cluster_cols=FALSE,
                        annotation_row =  ann,annotation_colors=ann_col,
                         clustering_method = "ward.D2",
                         color=mc,
                        border_color = "black")

pdf("DeepHB_selmotif_hm.pdf", width = 12, height = 16, pointsize = 10)
print(hp)
dev.off()


q()
## not used -------------
hp<- pheatmap::pheatmap(pat_mat,
                        clustering_method = "ward.D2",
                        border_color = NA)
motif_ord <- hp$tree_col
motif_ord <- motif_ord$labels[motif_ord$order]

seq_id_sel1 <- colnames(pat_mat)

fdr_df3 <- as.matrix(fdr_df2[,seq_id_sel1])
fdr_df3[fdr_df3>10] <- 10
hp<- pheatmap::pheatmap(fdr_df3,
                        clustering_method = "ward.D2",
                        border_color = NA)

pval_df <- reshape2::melt(fdr_df3)
head(pval_df)
colnames(pval_df) <- c("MP","Motifs","FDR")

mdf <- reshape2::melt(pat_mat)
head(mdf)
colnames(mdf) <- c("MP","Motifs","LSC")

row.names(mdf) <- paste0(mdf$MP,"_",mdf$Motifs)
row.names(pval_df) <- paste0(pval_df$MP,"_",pval_df$Motifs)
table(row.names(mdf) %in% row.names(pval_df))

mdf$FDR <- pval_df[row.names(mdf),"FDR"]
head(mdf)
rm(pval_df)

mdf$Motifs <- factor(mdf$Motifs, levels = rev(motif_ord))
mdf$MP <- factor(mdf$MP, levels = rev(mp_order))
mdf[mdf$FDR<3,"FDR"] <- 0

mdf$LSC <- pmax(mdf$LSC,0) #seqcounts
mdf[mdf$LSC<1,"FDR"] <- 0

p <- ggplot(mdf, aes(x=Motifs,y=MP, size=FDR))+
  geom_point(aes(color=LSC))+
  scale_size("negLog10FDR", range = c(0,10),breaks =c(0,5,10)) +
  #scale_color_gradient2(low = "blue",mid = "white", high = "red4",midpoint = 0)+
  scale_color_viridis_c(direction = -1,option = "F" ,
                        name = "Avg,\nAUC")+
  theme_minimal()+
  theme(#axis.line = element_line(colour = 'black',linewidth = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = 'black',linewidth=0.5),
    axis.text.x = element_text(colour = "black",size=5,
                               angle = 45,hjust=1),
    axis.text.y = element_text(face="bold",colour = "black",size=8),
    axis.title=element_blank())
png("DeepHB_motifs_MPenr_LSC.png", unit="in",width = 26, height = 12, res = 300)
print(p)
dev.off()
### homer #########
## obtain FDR matrix from homer run
fdr_df <- c()

for(x in mp_order){
  motif <- read.delim(paste0("TopicMP_CREs_57MP/motifenr_DH2/",x,"/knownResults.txt"),row.names = 1)
  motif <- motif[seq_ids,]
  del <- p.adjust(motif$P.value,method = "fdr")
  ##replace 0 fdr values with second minimum
  del2 <- min(del)
  if(del2==0){
    del3 <- min(del[del!=0])
    del[del==0] <- del3
    rm(del3)
  }
  
  fdr_df <- rbind(fdr_df,del)
  rm(motif,del,del2)
}
table(fdr_df==0)

fdr_df <- data.frame(fdr_df)
colnames(fdr_df) <- seq_ids
row.names(fdr_df) <- mp_order
fdr_df[1:4,1:4]

fdr_df2 <- (-1) *log10(fdr_df)
min(fdr_df2)
dim(fdr_df2)
fdr_df2[1:4,1:4]