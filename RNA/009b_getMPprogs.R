## Here individual LDA factors obtained from runs for all class+splits and ranks
## are processed and merged together to obtain the meta-gene program sets
suppressPackageStartupMessages({
  library(Matrix)
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

###---------
## 1. load samples-----
genelist <- list()
geneimplist <- list() ##need this for importance score later for merging
class_names <-readLines("RNA_Class.txt")
class_names <- sort(class_names)
lda_ranks <- seq(3,30,by=3)
class_names_wspl <-  c(outer(class_names, c("_A","_B"), paste0))
class_names_wspl <- sort(class_names_wspl) ##helps in annotation later

#class_names <- class_names[1:3]

##LDA--------
##loop on class names
for(idx in class_names_wspl){
  ##loop on lda ranks
  for(sel_rank in lda_ranks){
    id_name <- paste0(idx,":r",sel_rank)
    #load ldas data
    load(paste0(idx,"_",sel_rank,"_LDA.RData"))
    ldas <- colnames(rd_lda)
    for(ldaid in 1:length(ldas)){
      del <- data.frame(Score=sort(H[ldaid,], decreasing = TRUE)[1:100])## used 50, 100 and 150 for different runs
      del$Gene <- row.names(del)
      del$Rank <- seq(dim(del)[1])
      genelist[[paste0(id_name,"_",ldaid)]]  <- row.names(del)
      geneimplist[[paste0(id_name,"_",ldaid)]] <- del
      rm(del)
    }
    rm(ldaid,id_name,ldas,rd_lda,H,plot.data)
    
  }
  rm(sel_rank)
  
}
rm(idx)
#lengths(genelist)

## 3. similarity of MP across all MP----------
all_lda <- names(genelist)


# ##vectorized, same output as above
OS_AvB <- outer(all_lda, all_lda, FUN = Vectorize(function(i, j) {
  overlapS(genelist[[i]], genelist[[j]])
}))
rownames(OS_AvB) <- colnames(OS_AvB) <- all_lda


##vectorized jaccard
JS_AvB <- outer(all_lda, all_lda, FUN = Vectorize(function(i, j) {
  jaccard(genelist[[i]], genelist[[j]])
}))
rownames(JS_AvB) <- colnames(JS_AvB) <- all_lda

##ldas
save(JS_AvB, geneimplist, genelist, 
     class_names_wspl,class_names, file = "JS_AvBforLDAngene100.RData")
save(OS_AvB, geneimplist, genelist, 
     class_names_wspl,class_names, file = "OS_AvBforLDAngene100.RData")
q()
# post similarity matrix-----------

## Removing redundant lda and classnsplit specific lda that could represent arbitrary programs
## Basically, i first find ldas that share overlap with lda of other classnsplit
## The reason to split the class is to obtain class specific lda that are still robust and not random selection of ldas
## I then sorted them by overlap, the more overlap the more I need to select it.
## Then I remove other lda of the same classnsplit that also overlap with this selected lda
## If a lda in classnsplit doesn't overlap with other classnsplit I remove it too, this is to remove arbitrary lda,
## that don't represent robust programs
## This way i identified lda progroms that are shared across sample and identify a shared meta-gene program
## The chosen cut-off of 20%, 50% and then 33% are bit subjective
## the first 20% is to make sure these topcis overap with interesting feature across other classnsplits
## the second 50% is to make sure this lda atleast appear twice in another rank
## the third 33% is to remove redundant programs across rnaks
## from here -----------
load("OS_AvBforLDAngene100.RData")
lda_ranks <- seq(3,30,by=3)
##165 by 2 ldas per class
ann <- data.frame(lda=names(genelist),
                  Dataset=rep(class_names,each=330),
                  Rank=rep(lda_ranks,lda_ranks))
row.names(ann) <- ann$lda
ann$lda <- NULL

## the aim is to identify robust and non-redundant lda per class and obtain
## a union of this list

to_keep <- c()

for(sel_class in class_names){
  
  ## first the aim is robust lda per rank
  ## for this per rank, I will obtain lda that represented outside that rank
  ## overlap quite highly with another lda in the same rank, and then remove
  ## redudancy within a rank
  
  ## obtain lda associated with a class across both splits
  all_lda <- row.names(ann)[ann$Dataset==sel_class]
  ## subset the OS matrix to class
  OS_AvB_CT <- OS_AvB[all_lda,all_lda]
  dim(OS_AvB_CT)
  ## obtain all the ranks
  ranks <- paste0(":r",lda_ranks,"_") 
  
  ## this is to obtain robust non-redundant lda for each rank across both splits
  to_keep_a <- c() ## list of robust rank specific lda
  for(x in ranks){
    ## all lda for a rank across both split
    sel_lda <- grep(x,all_lda)
    ##extract OS of a rank with all the other ranks
    delOS <- OS_AvB_CT[sel_lda,-c(sel_lda)]
    ## max similarity for each ranks lda with other lda, sorted by max overlap
    delOS2 <- sort(apply(delOS, 1, max), decreasing = TRUE)
    ##only keep robust lda for rank, 
    ## these lda overlap with at least one another lda from other rank, CO=20%
    delOS2 <- delOS2[delOS2>=0.2] #.2 for OS, .11 for JS
    sely <- names(delOS2) ##robust part 1
    
     ## remove redundant lda within the same rank
    delOS3 <- OS_AvB_CT[sely,sely]
    delOS4 <- matrix(0, nrow = dim(delOS3)[1],ncol = dim(delOS3)[1])
    row.names(delOS4) <- colnames(delOS4) <- row.names(delOS3)
    ##create a matrix where overlaping lda are marked by 1
    delOS4[delOS3>=0.5] <- 1 ##0.5 for OS 0.333 J
    ## here the diag should not be set to 0, as it will remvoed in code below
    
    ##now removing  lda from a rank that are overlapped by selected
    ## robust lda ranked by overlap across other ranks
    sel_lday <- c()
    while (length(sely)>0) {
      if(length(sely)>1){
        ## take the first one
        selz <- sely[1]
        ## add it to selected lda
        sel_lday <- c(sel_lday,selz)
        ##extract OS for other lda which overlap with selected lda
        delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
        ##remove other lda from that overlap the selected lda from the selx
        delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
        sely <- sely[-c(which(sely %in% delOS5))]
        #print("reminaing lda :")
        #print(selx)
      }else {
        sel_lday <- c(sel_lday,sely)
        sely <- c()
      } 
    }
    
    
    to_keep_a <- c(to_keep_a,sel_lday)
    rm(sel_lda,delOS,delOS2,selx,sely,delOS3,delOS4,sel_lday,delOS5)
  }
  rm(x)
  
  ##now for each class remove the redundancy across ranks
  
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
  ##create a matrix where overlaping lda are marked by 1
  delOS4[delOS3>=0.5] <- 1 ##0.3 for OS ,0.215 for JS
  
  ##now removing  lda from a class that are overlapped by selected
  ## robust lda ranked by overlap across other ranks
  sel_lday <- c()
  while (length(sely)>0) {
    if(length(sely)>1){
      ## take the first one
      selz <- sely[1]
      ## add it to selected lda
      sel_lday <- c(sel_lday,selz)
      ##extract OS for other lda which overlap with selected lda
      delOS5 <- colnames(delOS4)[delOS4[selz,]==1]
      ##remove other lda from that overlap the selected lda from the selx
      delOS4 <- delOS4[-c(which(row.names(delOS4) %in% delOS5)),-c(which(row.names(delOS4) %in% delOS5))]
      sely <- sely[-c(which(sely %in% delOS5))]
      #print("reminaing lda :")
      #print(selx)
    }else {
      sel_lday <- c(sel_lday,sely)
      sely <- c()
    } 
  }
  
  to_keep <- c(to_keep,sel_lday)
  rm(all_lda, OS_AvB_CT, OS_AvB_CTx,ranks,to_keep_a)
  rm(row_max,sely,delOS3,delOS4,sel_lday,delOS5)
}
rm(sel_class)

length(to_keep) ## selected robust gene-sets

### now using jaccard matrix of selected gene-set for clustering

load("JS_AvBforLDAngene100.RData")
OS <- JS_AvB[to_keep,to_keep]
hp2 <-pheatmap::pheatmap(OS,
               clustering_method = 'complete',
               show_rownames = FALSE,show_colnames = FALSE)

hp2 <-pheatmap::pheatmap(OS,
               clustering_method = 'complete',
               #show_rownames = FALSE,
               show_colnames = FALSE,fontsize = 5)
del2 <- hp2$tree_col
del3 <- del2$labels
del3 <- del3[del2$order]
cluster = data.frame(X=cutree(hp2$tree_row, k = 40))
head(cluster)

ann$X <- "X"
annx <- ann[to_keep,]
annx$MP <- paste0("MP",cluster$X)
annx$X <- NULL
annx$Rank <- NULL
write.table(annx, file = "ann_40com_ldaMP100g_intial.txt", sep = "\t", quote = FALSE)

## no merging the obtained gene-set clusters to obtain meta-gene programs ---------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
})
load("OS_AvBforLDAngene100.RData")
annx <- read.delim("ann_40com_ldaMP100g_intial.txt",row.names = 1)
#annx <- read.delim("lda100g_mpfix.txt",row.names = 1) ## this is when I merge programs iteratively

mp_genes <- list()
ldampclust <- unique(annx$NewMP)

for(mpx in ldampclust){
  #identify lda clustered together
  sel_lda <- row.names(annx)[annx$NewMP==mpx]
  ##subset the list
  geneimplist_sub <- list()
  for(ldax in sel_lda){
    geneimplist_sub[[ldax]] <- geneimplist[[ldax]]
    
  }
  rm(ldax)
  ##merge the obtaind list of data frames
  merged_gene_df <- bind_rows(geneimplist_sub) %>%
    group_by(Gene) %>%
    summarise(Rank = mean(Rank, na.rm = TRUE), 
              count = n(), .groups = "drop")
  ##first rank by mean importance of genes and then by numebr of time gene appears
  ##this will result in most frequent gene and between gene that has same frequency 
  ##the one with higher mean importance/ or rank
  # merged_gene_df_ranked <- merged_gene_df %>%
  #   arrange(desc(count), desc(Score)) ##both descending
  merged_gene_df_ranked <- merged_gene_df %>%
    arrange(desc(count), Rank)
  mp_genes[[mpx]] <- merged_gene_df_ranked$Gene[1:100]
  
  rm(merged_gene_df,merged_gene_df_ranked,geneimplist_sub,sel_lda)
}
rm(mpx)
length(mp_genes)
lengths(mp_genes)


## checking overlap of the identified MP to further merging
all_mps <- names(mp_genes)
OS_AvB <- outer(all_mps, all_mps, FUN = Vectorize(function(i, j) {
  overlapS(mp_genes[[i]], mp_genes[[j]])
}))
rownames(OS_AvB) <- colnames(OS_AvB) <- all_mps

OS2 <- matrix(0,length(all_mps),length(all_mps))
row.names(OS2) <- colnames(OS2) <- all_mps
OS2[OS_AvB>0.3] <- 1
## heatmap
mb2 <- ((seq(1:101)-1)/100)
mc3 <- colorRampPalette(c("white","grey87","black"))(100)

hp2 <-pheatmap::pheatmap(OS2,
                         #clustering_method = 'ward.D2',
                         cluster_rows = FALSE,cluster_cols = FALSE,
                         breaks = mb2, color = mc3, fontsize = 8)
pdf("LDA_Orig_OS_30.pdf", width = 8, height = 8, pointsize = 10)
print(hp2)
dev.off()


##save this list of TopicMP CREs
save(mp_genes, file = "LDAMPgene_18MPfix_100g.RData") ## final set of 18 meta-gene programs

q()