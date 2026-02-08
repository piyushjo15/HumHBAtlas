### finding average expression of HOX gene targets in GC vs PN
suppressPackageStartupMessages({
  library(ggsci)
  library(ggrepel)
  library(tidyverse)
  library(pheatmap)
  library(Matrix)
})
DIR <- "~/PNana/"
setwd(DIR)
## ---
load("~/RNA/lgcountsHB.RData")
load("~/RNA/plotdataHBv2.RData")


## not subsetting

cl_lv <- c("Neuroblasts_I:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
plot.data <- plot.data[plot.data$Cluster %in% cl_lv,]

mdt <- mdt[,row.names(plot.data)]
dim(mdt)

cand <- c("TCF7L1","TCF7L2","ELF1","HMG20A",
          "SREBF2","HOXA3","HOXB2","HOXB3","HOXB4","HOXD3","PAX6")
plot.data$KNN_cl <- plot.data$Cluster



## hox genes
load("~/SCENIC/PN/metadata/PN_gmt.RData")
hox <- grep("HOX",names(regs),value = TRUE)
GSEA <- list()
for(x in hox){
  GSEA[[x]] <-regs[[x]]
}


### now for each of the hox genes I want to find the mean expression value of each targets
colx <- c("#e6f799","#BDDC32","#B2E688","#78AB51")
names(colx) <- c("Glioblasts:Progenitor_AES","Neuroblasts_I:A1_AES","GCUBC:GCP_PrN","GCUBC:GCP_PN")
colx <- colx[cl_lv]
expmat <- list()
for(x in hox){
  mat <- c()
  cand <- GSEA[[x]]
  for(y in cl_lv){
    cellsy <- row.names(plot.data)[plot.data$Cluster==y]
    mdtx <- mdt[cand,cellsy]
    del <- apply(mdtx,1,mean) ## mean is more real
    mat <- cbind(mat,del)
    rm(del)
  }
  colnames(mat) <- cl_lv
  expmat[[x]] <- mat
  rm(mat,cand)
  
}
pd <- c()
for(x in hox){
  del <- expmat[[x]]
  del2 <- reshape2::melt(del)
  del2$Gene <- x
  pd <- rbind(pd,del2)
  rm(del,del2)
}
head(pd)
pd$Var2 <- factor(pd$Var2, levels = rev(cl_lv))
p <- ggplot(pd, aes(y=Var2,x=value,fill=Var2))+
geom_boxplot(linewidth = 0.3,     # thinner outlines & whiskers
             fatten = 1.2,outlier.size = 1,
             outlier.alpha = 0.8)+
  theme_classic()+
  scale_fill_manual(values = colx)+
  theme(axis.line = element_line(colour = 'black',linewidth=.5),
        axis.ticks = element_line(colour = 'black',linewidth=.5),
        axis.text.y = element_text(colour = "black",size=10),
        axis.text.x = element_text(colour = "black",size=10),
        axis.title=element_blank(),
        strip.background = element_blank(),     # remove the box
        strip.placement  = "outside",           # keep strip outside panels
        strip.text.y.left = element_text(angle = 90, size=8),
        legend.position = "none")+
  scale_x_continuous(position = "top") +
  facet_wrap(~Gene,nrow=1,scales = "free_x",, strip.position = "top")
pdf("HOXGeneMedian.pdf",width = 8,height = 2,pointsize = 10)
print(p)
dev.off()
## bubbple plot using module score-------------
### findinding modscore
load("pdPNGCUBC_GRNs_MS.RData")
load("pdPNGCUBC_MBNMF_MS.RData")
load("pdPNGC_DEGt100PNGC_MS.RData")
load("pdPNGC_RandPNGC_MS.RData")
cl_lv <- c("Neuroblasts_A:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_EM","GCUBC:GCP_PN")
plot.data <- pd[pd$Cluster%in%cl_lv,]
grns <- colnames(scr)
hox <- grep("HOX",grns,value = TRUE)

hox <- c("TCF7L1","TCF7L2","POU3F2","SMAD3","SMAD2")
hox <- paste0("PN:",hox)

## didn't not down sample data
scr <- scr[row.names(plot.data),hox]

scr <- scr[row.names(plot.data),]
## 2.1 signficance, making bubble plot based on AUC and FDR

## using pair-wise wilcox and finding average values
marks <- pairwiseWilcox(t(scr), groups=plot.data$Cluster,direction="up")
sets <- data.frame(A=marks$pairs$first,B=marks$pairs$second)
dd <- marks$statistics
names(dd) <- paste0(sets$A," vs ",sets$B)

pval_df <- c()

for(x in cl_lv){
  b <- sets[sets$A==x,"B"]
  df_del <- c()
  for(y in b){
    z <- paste0(x," vs ",y)
    del <- data.frame(dd[[z]])
    del$B <-y
    del$GRN <- row.names(del)
    row.names(del) <- NULL
    df_del <- rbind(df_del,del[,c("GRN","B","AUC","p.value")])
    rm(del,z)
  }
  rm(y)
  
  #before finding fisher method adjusted p-value, fix the 0 pvalue
  ## either to the next minimum value or 1e-100, which ever is smaller
  min_p <- min(df_del$p.value[df_del$p.value!=0])
  if(min_p<1e-100){
    df_del[df_del$p.value==0,"p.value"] <- min_p
  }else{
    df_del[df_del$p.value==0,"p.value"] <- 1e-100
  }
  df_del$Cluster <- x
  pval_df <- rbind(pval_df,df_del)
  rm(min_p,df_del,b)
}
rm(x)
head(pval_df)

## this is to combine pvalue from various pair-wise test
combine_fisher <- function(p) {
  stat <- -2 * sum(log(p))
  pchisq(stat, df = 2*length(p), lower.tail = FALSE)
}

alpha <- 1e-5
pd <- pval_df %>%
  group_by(GRN, Cluster) %>%
  summarise(
    mean_AUC = mean(AUC, na.rm = TRUE),## mean AUC of a GRN in a class
    prop_sig = mean((AUC > 0.51) & (p.adjust(p.value, method = "bonferroni") < alpha)),## proportion of significant enrcihments
    fisher_p   = combine_fisher(p.value),##combined pvalue for this GRN for this Class
    Sig = ifelse(prop_sig > .5,"Yes","No"),## is it signficantly enriched?
    .groups = "drop"
  )

pd <- data.frame(pd)
head(pd)
pd[pd$mean_AUC<.1,"mean_AUC"] <- .1
pd$FDR <- p.adjust(pd$fisher_p,method = "bonferroni") # final correction of fisher_p
pd[pd$FDR==0,"FDR"] <- min(pd$FDR[pd$FDR!=0])/10 ## take care of zeros
#pd[pd$FDR<1e-10,"FDR"] <- 1e-10
pd$FDR <- -log10(pd$FDR)
pd[pd$mean_AUC<0.5,"Sig"] <- "No" ## if mean_AUC was <0.5 then remove from signficance
pd$Sig <- factor(pd$Sig,levels=c("No","Yes"))

pd$GRN <- factor(pd$GRN, levels =sort(hox))
pd$GRN <- factor(pd$GRN, levels =hox)

pd$Cluster <- factor(pd$Cluster, levels = cl_lv)

p<- ggplot(pd, aes(y = GRN, x = Cluster,color=Sig,
                   size = mean_AUC, fill = mean_AUC)) +
  geom_point(shape = 21) +
  scale_color_manual(values = c("grey87","black"))+
  #scale_size(range = c(2, 12),  name = "Log10FDR") +
  scale_size(range = c(1, 12),limits = c(0.1,1),  name = "Mean AUC") +
  #scale_fill_viridis_c(option = "blues",direction = -1, name = "Mean AUC") +
  scale_fill_gradientn(colors = c("white", "lightblue", "turquoise4"),    
                       limits = c(0.1, 1),
                       name = "Mean AUC") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(color="black",size=14),
    panel.grid.major = element_line(color="grey50",linewidth = .5),
    panel.grid.minor = element_blank()
  )+
  scale_x_discrete(position = "top")

w=.27*length(cl_lv)+4
h=.2*length(hox)+3

png(paste0("HOX_MS_bub_Wilcox_pw.png"),units = "in",width = w,height = h,res = 300)
print(p)
dev.off()

pdf(paste0("HOX_MS_bub_Wilcox_pw.pdf"),width = w,height = h,pointsize = 10)
print(p)
dev.off()


pdf(paste0("NonHOX_MS_bub_Wilcox_pw.pdf"),width = w,height = h,pointsize = 10)
print(p)
dev.off()

pdf(paste0("DEGt100_MS_bub_Wilcox_pw.pdf"),width = w,height = h,pointsize = 10)
print(p)
dev.off()

pdf("Rand100_MS_bub_Wilcox_pw.pdf",width = w,height = h,pointsize = 10)
print(p)
dev.off()
rm(pd,marks)


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

for(x in cl_lv){
  keep1 <- plot.data$Cluster==x
  msx <- ms2x[keep1,]
  del <- colMeans(msx)
  msd2 <- rbind(msd2,del)
  rm(keep1,msx,del)
}
rm(x)
rownames(msd2) <- cl_lv
msd2[1:4,1:4]


myColor <- colorRampPalette(c("deepskyblue4","grey93","magenta3"))(101)
# myBreaks <- c(seq(min(msd2), 0, length.out=ceiling(100/2) + 1),
#               seq(max(msd2)/100, max(msd2), length.out=floor(100/2)))
myBreaks <- c(seq(-2, 0, length.out=ceiling(100/2) + 1),
              seq(2/100, 2, length.out=floor(100/2)))
hp <- pheatmap(t(msd2), 
               clustering_method = "ward.D2",
               #clustering_distance_cols = "correlation",
               #clustering_distance_rows = "correlation",
               #annotation_row = cluster,
               annotation_colors = ann_col,
               #cluster_cols = FALSE,
               color=myColor, breaks=myBreaks,
               border_color = NA,
               fontsize=5)


marks <- findMarkers(t(scr),groups=plot.data$Cluster,direction="up",test.type="wilcox",pval.type="some")
marks$`Neuroblasts_A:A1_AES`
marks$`Glioblasts:Progenitor_AES`

###
ms2x <- data.frame(scr)
ms2x$Lv2 <- plot.data$Cluster
pd <- reshape2::melt(ms2x)
head(pd)
colx <- c("#e6f799","#BDDC32","#B2E688","#78AB51")
names(colx) <- c("Glioblasts:Progenitor_AES","Neuroblasts_A:A1_AES","GCUBC:GCP_EM","GCUBC:GCP_PN")
colx <- colx[cl_lv]

pd$Lv2 <- factor(pd$Lv2, levels = rev(cl_lv))
p <- ggplot(pd, aes(x=Lv2,y=value,fill=variable))+
  geom_boxplot(linewidth = 0.3,     # thinner outlines & whiskers
               fatten = 1.2,outlier.size = 1,
               outlier.alpha = 0.8)+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth=.5),
        axis.ticks = element_line(colour = 'black',linewidth=.5),
        axis.text.y = element_text(colour = "black",size=10),
        axis.text.x = element_text(colour = "black",size=10),
        axis.title=element_blank(),
        strip.background = element_blank(),     # remove the box
        strip.placement  = "outside",           # keep strip outside panels
        strip.text.y.left = element_text(angle = 90, size=8))

###

