## This scirpt obtained meta-regulatory program peak-set from the topic peaks
## obtained per class per splt across ranks

suppressPackageStartupMessages({
  library(pheatmap)
  library(Matrix)
  library(RColorBrewer)
})

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (round(intersection/union,4))
}
overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (round(intersection/minlen,4))
}
DIR <- "~/ATACana/Outs/TopicMP/"
setwd(DIR)
########

## 1. load samples-----
crelist <- list()
crenimplist <- list() ##need this for importance score later for merging
class_names <- readLines("~/RNA/RNA_Class.txt")
class_names <- sort(class_names)
topic_ranks <- seq(10,75,by=5)
class_names_wspl <-  c(outer(class_names, c("_A","_B"), paste0))
class_names_wspl <- sort(class_names_wspl) ##helps in annotation later
#class_names <- class_names[1:3]

## 2. use for loop to load topics for each class and each split ---------
##loop on class names
for(idx in class_names_wspl){
  ##loop on topic ranks
  for(sel_rank in topic_ranks){
    id_name <- paste0(idx,":r",sel_rank)
    #read topics
    del_topic <- read.delim(paste0(idx,"/topics/topic",sel_rank,"_regions.txt"))
    ntopics <- unique(del_topic$DataFrame_Name)
    for(topicid in ntopics){
      crelist[[paste0(id_name,"_",topicid)]] <- del_topic[del_topic$DataFrame_Name==topicid,"Row_Index"]
      del2 <- del_topic[del_topic$DataFrame_Name==topicid,c("Row_Index","Value")]
      del2$Rank <- seq(dim(del2)[1])
      row.names(del2) <- del2$Row_Index
      crenimplist[[paste0(id_name,"_",topicid)]] <- del2
      rm(del2)
    }
    rm(topicid,del_topic,id_name,ntopics)
    
  }
  rm(sel_rank)
  
}
rm(idx)

## 3. similarity of MP across all MP----------
all_topics <- names(crelist)


## vectorized upper triangle
n <- m <- length(all_topics)
score_mat <- outer(
  1:n, 1:m,
  Vectorize(function(i, j) {
    if (i <= j) overlapS(crelist[[i]], crelist[[j]]) else NA
  })
)
OS_AvB <- score_mat
OS_AvB[lower.tri(OS_AvB)] <- t(score_mat)[lower.tri(score_mat)]
rm(score_mat)
rownames(OS_AvB) <- colnames(OS_AvB) <- all_topics
save(OS_AvB, crenimplist, crelist, class_names_wspl,class_names, file = "OS_AvBfortopic.RData")


score_mat <- outer(
  1:n, 1:m,
  Vectorize(function(i, j) {
    if (i <= j) jaccard(crelist[[i]], crelist[[j]]) else NA
  })
)
JS_AvB <- score_mat
JS_AvB[lower.tri(JS_AvB)] <- t(score_mat)[lower.tri(score_mat)]
rm(score_mat)
rownames(JS_AvB) <- colnames(JS_AvB) <- all_topics
#JS_AvB[1:4,1:4]
#OS_AvB[1:4,1:6]
save(JS_AvB, crenimplist, crelist, class_names_wspl,class_names, file = "JS_AvBfortopic.RData")
q()

## 4. Redundancy ----------
## removing redundant topics and classnsplit specific topics that could represent arbitrary programs
## Basically, i first find Topics that share overlap with Topics of other classnsplit
## The reason to split the class is to obtain class specific topics that are still robust and not random selection of cres
## I then sorted them by overlap, the more overlap the more I need to select it.
## Then I remove other Topics of the same classnsplit that also overlap with this selected Topic.
## If a Topic in classnsplit doesn't overlap with other classnsplit I remove it too, this is to remove arbitrary topics,
## that don't represent robust programs
## This way i identified Topic progroms that are shared across sample and identify a shared meta Topicprogram
## The chosen cut-off of 20%, 50% and then 33% are bit subjective
## the first 20% is to make sure these topcis overap with interesting feature across other classnsplits
## the second 60% is to make sure this topic atleast appear twice in another rank
## the third 33% is to remove redundant programs but since these are CREs and not genes
## I am using a bit higher overlap cut-off for redundancy
### form here ----------
load("OS_AvBfortopic.RData")
topic_ranks <- seq(10,75,by=5)
ann <- data.frame(topics=names(crelist),
                  Dataset=rep(class_names,each=1190),
                  Rank=rep(topic_ranks,topic_ranks))
## reduce number of topics
## this is for v2
topic_ranks <- seq(10,50,by=5)
ann1 <- ann
ann <- ann1[ann1$Rank %in% topic_ranks,]

#
row.names(ann) <- ann$topics
ann$topics <- NULL
colclass <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
colx <- colclass$Color
names(colx) <- row.names(colclass)
rm(crenimplist,crelist)

## the aim is to identify robust and non-redundant topics per class and obtain
## a union of this list

to_keep <- c()

for(sel_class in class_names){
  
  ## first the aim is robust topics per rank
  ## for this, per rank, I will obtain topics that represented outside that rank
  ## overlap quite highly with another topic in the same rank, and then remove
  ## redudancy within a rank
  
  ## obtain topics associated with a class across both splits
  all_topics <- row.names(ann)[ann$Dataset==sel_class]
  ## subset the OS matrix to class
  OS_AvB_CT <- OS_AvB[all_topics,all_topics]
  dim(OS_AvB_CT)
  ## obtain all the ranks
  ranks <- paste0(":r",topic_ranks,"_") 
  
  ## this is to obtain robust non-redundant topics for each rank across both splits
  to_keep_a <- c() ## list of robust rank specific topics
  for(x in ranks){
    ## all topics for a rank across both split
    sel_topics <- grep(x,all_topics)
    ##extract OS of a rank with all the other ranks of same class
    delOS <- OS_AvB_CT[sel_topics,-c(sel_topics)]
    ## max similarity for each ranks topics with other topics, sorted by max overlap
    delOS2 <- sort(apply(delOS, 1, max), decreasing = TRUE)
    ##only keep robust topics for rank, 
    ## these topics overlap with at least one another nmf from other rank, CO=20%
    delOS2 <- delOS2[delOS2>=0.2] #.2 for OS, .11 for JS
    selx <- names(delOS2) ##robust part 1
    
    ## This is needed for ATAC as I have more ranks so there should be robust
    ## topics within same rank across splits
    ## now for these ranks specific robust topics, further identify robust topics
    # that are 33% overlapped by atleast one more topics from the same rank, maybe
    # in the other split. I also remvoe these topics to remove redundancy
    ## now extract OS of robust topics for a rank
    delOS3 <- OS_AvB_CT[selx,selx]
    delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
    row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
    ##create a matrix where overlaping topics are marked by 1
    delOS4[delOS3>=0.33] <- 1 #.6 for OS, 0.43 for JS
    diag(delOS4) <- 0 # don't self-count
    rs <- rowSums(delOS4)
    keep <- rs>0
    sely <- row.names(delOS3)[keep] ##robust part 2
    
    ## remove redundant topics within the same rank
    ##sely <- selx ### skipping part 2
    delOS3 <- OS_AvB_CT[sely,sely]
    delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
    row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
    ##create a matrix where overlaping topics are marked by 1
    delOS4[delOS3>=0.33] <- 1 ##0.5 for OS 0.333 J
    ## here the diag should not be set to 0, as it will remvoed in code below
    
    ##now removing  topics from a rank that are overlapped by selected
    ## robust nmf ranked by overlap across other ranks
    sel_topicsy <- c()
    while (length(sely)>0) {
      if(length(sely)>1){
        ## take the first one
        selz <- sely[1]
        ## add it to selected nmf
        sel_topicsy <- c(sel_topicsy,selz)
        ##extract OS for other topics which overlap with selected nmf
        delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
        ##remove other topics from that overlap the selected nmf from the selx
        delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
        sely <- sely[-c(which(sely %in% delOS5))]
        #print("reminaing topics :")
        #print(selx)
        rm(selz)
      }else {
        sel_topicsy <- c(sel_topicsy,sely)
        sely <- c()
      } 
    }
    
    
    to_keep_a <- c(to_keep_a,sel_topicsy)
    rm(sel_topics,delOS,delOS2,selx,sely,delOS3,delOS4,sel_topicsy,delOS5)
  }
  rm(x)
  
  ##now for each class remove the redundancy across ranks
  
  ## here I am ranking them based on overlapp across clasess
  ## but this can select again uniquesness in a class
  #keep2 <-row.names(ann)[ann$Dataset!=sel_class]
  #OS_AvB_CTx <- OS_AvB[to_keep_a,keep2]
  
  ## this is ranked by robust in a class
  ## this looks better as the final heamtap looks more class 
  ## specific
  OS_AvB_CTx <- OS_AvB[to_keep_a,to_keep_a]
  
  diag(OS_AvB_CTx) <-0
  row_max <- apply(OS_AvB_CTx,1,max)
  row_max <- sort(row_max,decreasing = TRUE) ## these are also sorted in ascending rank
  sely <- names(row_max)
  
  ##repeating above method
  delOS3 <- OS_AvB_CT[sely,sely]
  delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
  row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
  ##create a matrix where overlaping topics are marked by 1
  delOS4[delOS3>=0.66] <- 1 ##0.3 for OS ,0.215 for JS
  
  ##now removing  topics from a class that are overlapped by selected
  ## robust topic ranked by overlap across other ranks
  sel_topicsy <- c()
  while (length(sely)>0) {
    if(length(sely)>1){
      ## take the first one
      selz <- sely[1]
      ## add it to selected nmf
      sel_topicsy <- c(sel_topicsy,selz)
      ##extract OS for other topics which overlap with selected nmf
      delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
      ##remove other topics from that overlap the selected nmf from the selx
      delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
      sely <- sely[-c(which(sely %in% delOS5))]
      #print("reminaing topics :")
      #print(selx)
    }else {
      sel_topicsy <- c(sel_topicsy,sely)
      sely <- c()
    } 
  }
  
  to_keep <- c(to_keep,sel_topicsy)
  rm(all_topics, OS_AvB_CT, OS_AvB_CTx,ranks,to_keep_a)
  rm(row_max,sely,delOS3,delOS4,sel_topicsy,delOS5)
}
rm(sel_class)

length(to_keep)


load("JS_AvBfortopic.RData")
rm(crenimplist,crelist)
#OS <- OS_AvB[to_keep,to_keep]
OS <- JS_AvB[to_keep,to_keep]


## 5. Identifying Meta Topic groups ------------
# pheatmap apprach --------
hp2 <-pheatmap::pheatmap(OS,
                         clustering_method = 'complete',
                         show_rownames = FALSE,show_colnames = FALSE)


## after obtaining the selected programs, now obtaining MP definition of their clusters
del2 <- hp2$tree_col
del3 <- del2$labels
del3 <- del3[del2$order]
cluster = data.frame(X=cutree(hp2$tree_row, k = 100)) ## 75 for v1, v2 and 100
head(cluster)

ann$X <- "X"
annx <- ann[to_keep,]
annx$MP <- paste0("MP",cluster$X)
annx <- annx[del3,]
annx$X <- NULL
annx$Rank <- NULL


mpcx <-sample(colorRampPalette(pal_igv()(51))(length(unique(annx$MP))))
names(mpcx) <- unique(annx$MP)

ann_col <- list(Dataset=colx,
                MP=mpcx)

## heatmap
mb <- ((seq(1:101)-1)/100)/2
mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)

hp2 <-pheatmap::pheatmap(OS[del3,del3],
                         #clustering_method = 'complete',
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         annotation_col = annx,
                         annotation_colors = ann_col,
                         #show_rownames = FALSE,
                         show_colnames = FALSE,
                         breaks = mb, color = mc2, fontsize = 5)
pdf("LDA_100com.pdf", width = 20, height = 20, pointsize = 10)
print(hp2)
dev.off()
write.table(annx, file = "ann_75com_ldaMP.txt", sep = "\t", quote = FALSE)


## since there are so many data points, I am subsetting the heatmap to see better

## fix, merge some clusters, unmerge others
annx <- read.delim("ann_75com_ldaMPfix.txt",row.names = 1)

nmpx <-sample(colorRampPalette(pal_igv()(51))(length(unique(annx$MP))))
names(nmpx) <- unique(annx$MP)
nmpx2 <-sample(colorRampPalette(pal_igv()(51))(length(unique(annx$NewMPv2))))
names(nmpx2) <- unique(annx$NewMPv2)
colclass <- read.delim("~/MBsnANA/Diemana/scepro/Classcol.txt",row.names = 1)
colx <- colclass$Color
names(colx) <- row.names(colclass)
ann_col <- list(Class=colx,
                NewMPv2=nmpx2,
                MP=nmpx)

## full fixed heatmap
mb <- ((seq(1:101)-1)/100)/2
mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)
anny <- annx
mpcl <- unique(annx$MPCL)
anny$NewMP <- anny$NewMP_old <- anny$MPCL <- NULL
for(x in mpcl){
  del <- row.names(annx)[annx$MPCL==x]
  hpx <-pheatmap::pheatmap(OS[del,del],
                           cluster_rows = FALSE,cluster_cols = FALSE,
                           annotation_col = anny,
                           annotation_colors = ann_col,
                           #show_rownames = FALSE,
                           show_colnames = FALSE,border_color = NA,
                           breaks = mb, color = mc2, fontsize = 5)
  pdf(paste0("LDA_75com_",x,".pdf"), width = 12, height = 12, pointsize = 10)
  print(hpx)
  dev.off()
  rm(del,hpx)
  
}
anny$Old_NewMPv2 <- anny$MP <- NULL

hp2 <-pheatmap::pheatmap(OS,
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         annotation_col = anny,
                         annotation_colors = ann_col,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         breaks = mb, color = mc2, fontsize = 5)
pdf("LDA_75comfix_nocluster.pdf", width = 20, height = 20, pointsize = 10)
print(hp2)
dev.off()
pdf("LDA_70comfix_finalOS.pdf", width = 20, height = 20, pointsize = 10)
print(hp2)
dev.off()
## 5. merging TOpic program to obtain top2k CREs per Topic MPs in a cluster-------
## In this step, before obtaining final CREs for each MP, I am doing analysis
## to identify which MPs can be merged

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggsci)
})
overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (round(intersection/minlen,4))
}
load("OS_AvBfortopic.RData")

annx <- read.delim("ann_75com_ldaMPfix.txt",row.names = 1)
to_keep <- row.names(annx)
OS <- OS_AvB[to_keep,to_keep]

mp_cres <- list()

topicmpclust <- unique(annx$NewMPv2)

for(mpx in topicmpclust){
  #identify topics clustered together
  sel_topics <- row.names(annx)[annx$NewMPv2==mpx]
  ##subset the list
  crenimplist_sub <- list()
  for(topicx in sel_topics){
    crenimplist_sub[[topicx]] <- crenimplist[[topicx]]
    
  }
  rm(topicx)
  ##merge the obtaind list of data frames
  ##merge the obtaind list of data frames
  merged_cre_df <- bind_rows(crenimplist_sub) %>%
    group_by(Row_Index) %>%
    summarise(Rank = mean(Rank, na.rm = TRUE), 
              count = n(), .groups = "drop")
  
  #select top 2k CREs ranked by importance
  merged_cre_df_ranked <- merged_cre_df %>%
    arrange(desc(count), Rank)
  ncres <- dim(merged_cre_df_ranked)[1]
  mp_cres[[mpx]] <- merged_cre_df_ranked$Row_Index[1:2500] ##2000 in v1
  
  rm(merged_cre_df,merged_cre_df_ranked,crenimplist_sub,sel_topics,ncres)
}
rm(mpx)
length(mp_cres)

##vectorized, same output as above
OS_AvB2 <- outer(topicmpclust, topicmpclust, FUN = Vectorize(function(i, j) {
  overlapS(mp_cres[[i]], mp_cres[[j]])
}))
rownames(OS_AvB2) <- colnames(OS_AvB2) <- topicmpclust

# ## heatmap
# mb2 <- ((seq(1:101)-1)/100)/2
# mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)
# 
# hp2 <-pheatmap::pheatmap(OS_AvB,
#                          cluster_rows = FALSE,cluster_cols = FALSE,
#                          #clustering_method = 'complete',
#                          #show_rownames = FALSE,
#                          show_colnames = FALSE,
#                          border_color = NA,
#                          breaks = mb2, color = mc2, fontsize = 5)
# pdf("LDA_75com_mergeNewMPv2.pdf", width = 10, height = 10, pointsize = 10)
# print(hp2)
# dev.off()
## overlap >=0.33 --> 1 else 0
## heatmap
OS2 <- matrix(0, length(topicmpclust),length(topicmpclust))
row.names(OS2) <- colnames(OS2) <- topicmpclust
OS2[OS_AvB2>=0.15] <- 1
mc2 <- colorRampPalette(c("white","black"))(100)

hp2 <-pheatmap::pheatmap(OS2,
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         #clustering_method = 'complete',
                         #show_rownames = FALSE,
                         show_colnames = FALSE,
                         border_color = NA,
                         color = mc2, fontsize = 5)

pdf("LDA_75com_mergeNewMP_15.pdf", width = 10, height = 10, pointsize = 10)
print(hp2)
dev.off()


head(sort(OS_AvB2["MP55",],decreasing = TRUE))
sort(OS["GC_B:r50_Topic44",],decreasing = TRUE)[5:30]

### finally obtaining list of cres per MP-----------------
## this is part two after further merging similar NewMPs
annx <- read.delim("ann_75com_ldaMPfix.txt",row.names = 1)
load("OS_AvBfortopic.RData")
mp_cres <- list()
topicmpclust <- unique(annx$NewMPv2)

for(mpx in topicmpclust){
  #identify topics clustered together
  sel_topics <- row.names(annx)[annx$NewMPv2==mpx]
  ##subset the list
  crenimplist_sub <- list()
  for(topicx in sel_topics){
    crenimplist_sub[[topicx]] <- crenimplist[[topicx]]
    
  }
  rm(topicx)
  ##merge the obtaind list of data frames
  ##merge the obtaind list of data frames
  merged_cre_df <- bind_rows(crenimplist_sub) %>%
    group_by(Row_Index) %>%
    summarise(Rank = mean(Rank, na.rm = TRUE), 
              count = n(), .groups = "drop")
  
  #select top 2k CREs ranked by importance
  merged_cre_df_ranked <- merged_cre_df %>%
    arrange(desc(count), Rank)
  ncres <- dim(merged_cre_df_ranked)[1]
  ## v2, 500 in ext
  if(ncres>=2500){
    mp_cres[[mpx]] <- merged_cre_df_ranked$Row_Index[1:2500]
  }else if(ncres>2000 & ncres <2500){
    mp_cres[[mpx]] <- merged_cre_df_ranked$Row_Index[1:ncres]
  }else {
    mp_cres[[mpx]] <- merged_cre_df_ranked$Row_Index[1:ncres]
  }
  
  rm(merged_cre_df,merged_cre_df_ranked,crenimplist_sub,sel_topics,ncres)
}

rm(mpx)
length(mp_cres)
length(mp_cres_ext)

lengths(mp_cres)
lengths(mp_cres_ext)

##save this list of TopicMP CREs
#save(mp_cres, mp_cres_ext, file = "TopiCMPCREslist_67MP.RData") ## old one used earlier now Old_NewMP2 annotation

save(mp_cres, file = "TopiCMPCREslist_57MP.RData") ## new on 040525
##newer version save each of them as bed file
## for each MP I am selecting 2500 CREs from top
## than out of 2500 CREs I am using random selection of 2k CRE for model
## and 500 CRE for an additional test.
## however, this shouldn't really make a huge difference as in the test portion
## the model never sees that either.
## need to think about this more.
load("TopiCMPCREslistv2.RData")
topicmpclust <- names(mp_cres)
topicmpclust2 <- names(mp_cres_ext)
sam <- sample(1:2000,500)

for(mpx in topicmpclust){
  del <- data.frame(CRE=mp_cres[[mpx]])
  del2 <- del %>% separate(CRE,c("A","B","C"))
  del2$D <- del$CRE
  del2$E <- mpx
  write.table(del2, paste0("TopicMP_CREs/",mpx,"_CREs.bed"),
              row.names=FALSE,
              col.names = FALSE, quote = FALSE, sep="\t")
  rm(del,del2)
  # if(mpx %in% topicmpclust2){
  #   cre <- mp_cres_ext[[mpx]]
  #   del <- data.frame(CRE=cre)
  #   del2 <- del %>% separate(CRE,c("A","B","C"))
  #   del2$D <- del$CRE
  #   del2$E <- mpx
  #   write.table(del2,paste0("TopicMP_CREs/",mpx,"_CREs.ext.v2.bed"),
  #               row.names=FALSE,
  #               col.names = FALSE, quote = FALSE, sep="\t")
  #   rm(del,del2,cre)
  # }
  
}
rm(mpx)

## new using extended 2500 CREs
load("TopiCMPCREslist_57MP.RData")
topicmpclust <- names(mp_cres)
for(mpx in topicmpclust){
  del <- data.frame(CRE=mp_cres[[mpx]])
  del2 <- del %>% separate(CRE,c("A","B","C"))
  del2$D <- del$CRE
  del2$E <- mpx
  
  ## sample rest 250 for test
  write.table(del2, paste0("TopicMP_CREs_57MP/",mpx,"_CREs.bed"),
              row.names=FALSE,
              col.names = FALSE, quote = FALSE, sep="\t")
  

  rm(del,del2)
}
rm(mpx)



## old --------------


#
# ##using hclust--------
# ## using OS matrix as cosine matrix
# dmat <- as.dist(1 - OS) # Using 1 - correlation as distance
# hc <- hclust(dmat, method = "ward.D2") ##this look good for now
# #plot(hc)
# 
# #hc
# cluster = data.frame(X=cutree(hc, k = 100))
# head(cluster)
# ann$X <- "X"
# annx <- ann[to_keep,]
# annx$MP <- paste0("MP",cluster$X)
# annx$X <- NULL
# mpx <-sample(colorRampPalette(pal_igv()(51))(length(unique(annx$MP))))
# names(mpx) <- unique(annx$MP)
# ann_col <- list(Dataset=colx,
#                 MP=mpx)
# 
# mb <- ((seq(1:101)-1)/100)/2
# mc2 <- colorRampPalette(c("white","sandybrown","brown4"))(100)
# library(circlize)
# col_fun = colorRamp2(c(0,0.25, 0.5), c("white", "sandybrown", "brown4"))
# col_fun(seq(0,1,by=0.1))
# 
# hp3 <- Heatmap(OS,
#                show_row_names = FALSE,
#                show_column_names = FALSE,
#                cluster_rows = hc,cluster_columns = hc,
#                col=col_fun)
# pdf("topics_50hclut.pdf", width = 16, height = 16, pointsize = 10)
# print(hp3)
# dev.off()
# write.table(annx, file = "ann_TopicMPhclustcom.txt", sep = "\t", quote = FALSE)
# annx <- read.delim("ann_TopicMPhclustcom.txt")
# 
# # tiff("Topicprotree_colv2.tiff", units="in", width=6, height=4, res=300)
# # plot(hp2$tree_col,cex=0.5 )
# # dev.off()

## old robust topic method --------
all_topics <- names(crelist)
##for Jaccards
OS_AvB <- JS_AvB
to_keep <- c()
for(x in class_names_wspl){
  ## all topics for a classnsplit
  sel_topics <- grep(x,all_topics)
  ##extract OS of a classnsplit with rest of classnsplit
  delOS <- OS_AvB[sel_topics,-c(sel_topics)]
  ## max similarity for each classnsplit topics with other topics, sorted by max overlap
  delOS2 <- sort(apply(delOS, 1, max), decreasing = TRUE)
  ##only keep robust topics for classnsplit, 
  ## these topics overlap with at least one another topic from other classnsplit, CO=20%
  delOS2 <- delOS2[delOS2>=0.2]
  selx <- names(delOS2)
  
  #now for these classnsplit specific robust topics, further identify robust topics
  # that are 50% overlapped by atleast one more topics from the same classnsplit 
  ## now extract OS of a classnsplit
  delOS3 <- OS_AvB[selx,selx]
  delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
  row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
  ##create a matrix where overlaping topics are marked by 1
  delOS4[delOS3>=0.5] <- 1
  diag(delOS4) <- 0
  rs <- rowSums(delOS4)
  keep <- rs>0
  
  #table(keep)
  ## further filter to remove redundant classnsplit specific topics
  ## ## and repeat above matrix
  sely <- row.names(delOS3)[keep]
  delOS3 <- OS_AvB[sely,sely]
  delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
  row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
  ##create a matrix where overlaping topics are marked by 1
  delOS4[delOS3>=0.33] <- 1
  
  ##now removing  topics from a classnsplit that are overlapped by selected
  ## robust topic ranked by overlap across other classnsplits
  sel_topicsy <- c()
  while (length(sely)>0) {
    if(length(sely)>1){
      ## take the first one
      selz <- sely[1]
      ## add it to selected topic
      sel_topicsy <- c(sel_topicsy,selz)
      ##extract OS for other topics which overlap with selected topic
      delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
      ##remove other topics from that overlap the selected topic from the selx
      delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
      sely <- sely[-c(which(sely %in% delOS5))]
      #print("reminaing topics :")
      #print(selx)
    }else {
      sel_topicsy <- c(sel_topicsy,sely)
      sely <- c()
    } 
  }
  
  to_keep <- c(to_keep,sel_topicsy)
  rm(sel_topics,delOS,delOS2,selx,sely,delOS3,delOS4,sel_topicsy,delOS5)
}
rm(x)
length(to_keep)






#

