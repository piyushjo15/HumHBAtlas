## This script merges overlapping grns from developmentally related super clusters

overlapS <- function(a, b) {
  intersection = length(intersect(a, b))
  minlen = min(length(a),length(b))
  return (intersection/minlen)
}

celltypes <- c("Oligopro","OligoMat","Astrocytes","Glioblasts")
sams <- c("OPC","MOG","AS","GB")
DIR = "~/MBsnANA/HBana/SCENICana/SCENICOut/Supercluster/"
setwd(DIR)
## 1. load grns -----------
df <- c()
GSEA <- list()
for(x in 1:length(celltypes)){
  y <- celltypes[x]
  del <- get(load(paste0(y,"/metadata/",y,"_gmt.RData")))
  a <- names(del)
  names(del) <- paste0(sams[x],"_",names(del))
  b <- names(del)
  
  GSEA <- c(GSEA,del)
  del2 <- data.frame(GRN=b,TF=a,len=lengths(del),CT=y)
  df <- rbind(df,del2)
  rm(y,del,a,b,del2)
}
rm(x)
grnx <- df$GRN[df$len<20]
for(x in grnx ){
  GSEA[[x]] <- NULL
}

df <- df[!(df$GRN %in% grnx),]
head(df)

## 2. find overlaps for TF-GRNS that share TFs--------
tfd <- df$TF[duplicated(df$TF)]
df2 <- df[df$TF %in% tfd,]

grns <- df2$GRN 
OS_AvB <- matrix(0, nrow = length(grns),ncol = length(grns))
row.names(OS_AvB) <- colnames(OS_AvB) <- grns

for(i in grns){
  for (j in grns){
    OS_AvB[i,j] <- overlapS(GSEA[[i]],GSEA[[j]])
  }
  rm(j)
}
rm(i)
OS_AvB[1:3,1:3]
mc2 <- colorRampPalette(c("white",RColorBrewer::brewer.pal(9,"Greys")))(100)

mb <- ((seq(1:101)-1)/100)

hp2 <-pheatmap::pheatmap(OS_AvB, clustering_method = 'ward.D2',
                         breaks = mb, color = mc2, fontsize = 8)

## 3. finding GRNS for same TFs that have less than half of novel member----
diag(OS_AvB) <- 0

##this method removes smaller GRNS that overlap with bigger GRNS
## but if there are three GRNS: A, B and C, in this order of size, and
## B overlaps A and C but A and C don't overlap much, it will only remove B and 
## not C, because after removing B, C is unique
torem <- c()
for(x in tfd){
  df_del <-df2[df2$TF==x,]
  ##sort GRNS by size
  df_del <- df_del[order(df_del$len,decreasing = TRUE),]
  OS_del <- OS_AvB[df_del$GRN,df_del$GRN]
  grns_del <- colnames(OS_del)
  ## remove non-unique GRN from set with smallest constituents
  for(y in grns_del){
    if(!(y %in% torem)){
      del1 <- grns_del[OS_del[y,]>=0.5]
      del1 <- del1[!(del1 %in% torem)] ##don't use GRNS already remvoed
      del2 <- c()
      if(length(del1)>0){
        for(z in del1){
          a <- df2[y,"len"]
          b <- df2[z,"len"]
          if(a>=b){
            del2 <- c(del2,z)
            
          } else {
            del2 <- c(del2,y)
          }
          rm(a,b)
        }
        rm(z)
      }
      torem <- c(torem,del2)
      rm(del1,del2)
    } 
    
  }
  rm(df_del,OS_del,grns_del,y)

}
rm(x)
torem <- unique(torem)

for(x in torem ){
  GSEA[[x]] <- NULL
}
df <- df[(df$GRN %in% names(GSEA)),]

save(GSEA, df, file = "MergedGlialGRNS.RData")

## for some plots---------
load("MergedGlialGRNS.RData")
df_del <- df[df$TF=="SOX4",]
sel <- df_del$GRN
GSEA2 <- list()
for(x in sel){
  GSEA2[[x]] <- GSEA[[x]]
}
names(GSEA2) <-  c("OPC" , "MOG " ,"Astrocytes","Glioblast")
library(ggvenn)
ggvenn(GSEA2, text_size = 8,fill_color=c("blue", "yellow", "green", "red"),
       show_percentage = FALSE)

