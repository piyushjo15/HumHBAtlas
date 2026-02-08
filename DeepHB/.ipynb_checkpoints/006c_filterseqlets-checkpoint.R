library(tidyverse)

setwd("~/DeepHB/Outs/")

## select patterns based on some filter
pat_names <- readLines("~/DeepHB/Outs/all_pat_names.txt")
seq_ids <- paste0("seqid_",sprintf("%03d", 1:165)) ## in total 165 patterns were obtained
pat_names <- data.frame(X=pat_names)
pat_names <- pat_names %>% separate(X,c("A","B"))
pat_names <- pat_names$B
pat_list <- str_split(pat_names,"")
keep <- lengths(pat_list)>3
table(keep)
seq_ids <- seq_ids[keep]

### identify imporant motifs from seqcount
## this is a rMP by motif/seqlet count matrix
load( "~/DeepHB/Outs/Seqletlogcounts_allPat.RData")
pattern_matrix <- pattern_matrix[,seq_ids]
mps <- row.names(pattern_matrix) ## these are 57 rMPs
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


## also identifying top pattern per MP

seq_id_list <- list()
seq_id_mat <- c()
seq_ids_sel <- c()
for(x in mps){
  del <- names(sort(pattern_matrix[x,],decreasing = TRUE))[1:2]
  seq_id_list[[x]] <- del
  seq_id_mat <- rbind(seq_id_mat,del)
  seq_ids_sel <- c(seq_ids_sel,del[1])
  rm(del)
}
rm(x)
row.names(seq_id_mat) <- mps
seq_ids_sel <- unique(seq_ids_sel) ## 65 in top 3, 43 in top 2, 28 in top 1
sel_seq_ids <- unique(unlist(seq_id_list))
length(sel_seq_ids)

row.names(seq_id_mat) <- mps
write.table(seq_id_mat, file = "top2seqidperMP.txt",sep="\t",quote = FALSE)

save(seq_id_list, file = "topseqidsperMP.RData")
write(seq_ids_sel, file = "topseqidperMP.txt")


write(sel_seq_ids, file ="~/DeepHB/Outs/selectedseqidsforplot.txt" )
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


q()
