## This script performs enrichment analysis of CRE sets
## in pontine nuclei or granule cell precursors
suppressPackageStartupMessages({
  library(ArchR)
  library(tidyverse)
  library(ComplexHeatmap)
  library(scater)
  library(scran)
  library(AUCell)
  
})

DIR_ATAC <- "~/ATACana/Outs/"
DIR_GRN <- "~/SCENIC/"

##
## GRN sets
load(paste0(DIR_GRN,"PN/metadata/PN_gmt.RData"))

hox <- grep("HOX",names(regs),value = TRUE)
non_hox <- c("TCF7L1","TCF7L2","POU3F2","SMAD3","SMAD2")
sel_cls <- c("Neuroblasts_I:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_PrN","GCUBC:GCP_PN")
########## ATAC part ###############
## part3. load the region sets-----------
## region sets
df <- read.delim(paste0(DIR_GRN,"PN/metadata/eRegulons_fix.txt"))
regCREs <- df[df$TF %in% hox,] ## or non_hox

#Finding HOX binding CREs
head(regCREs)
regCREs$X <- paste0(regCREs$TF,"_",regCREs$Region)
regCREs <- regCREs[!duplicated(regCREs$X),]
row.names(regCREs) <- regCREs$X

## cre-set for AUC
GSEA <- list()
for(x in hox){
  GSEA[[x]] <- regCREs[regCREs$TF==x,"Region"]
}
lengths(GSEA)

## Finding non-hox CREs of HOX targets per HOXgene
## For each HOX gene, find its targets, find all the CREs of those targets
## then remove HOX binding CREs

GSEA <- list()
GSEA_tf <- list()

for(x in hox){
  tars <- regCREs[regCREs$TF==x,"Gene"]
  del <- df[df$Gene %in% tars,]
  GSEA[[x]] <- unique(del[del$TF!=x,"Region"])
  GSEA_tf[[x]] <- unique(del[del$TF!=x,"TF"])
  
}
lengths(GSEA)
lengths(GSEA_tf)


## finding CREs associated with GC or PN genes based on SCENIC+ selected CREs
load("PNana/DEG_PNGCtop100.RData")
deg <- GSEA
rm(GSEA)
## only 71 GC and 78 Pontine genes have CREs
hox <- names(deg)
GSEA <- list()
for(x in hox){
  gx <- deg[[x]]
  GSEA[[x]] <- unique(df[df$Gene %in% gx,"Region"])
  
}
lengths(GSEA)

## part 4. ATAC part------------

## load plot metadata and CRE signature
load("~/ATACana/Outs/plotdataHBATACv6.RData")
plot.data <- plot.data[plot.data$Cluster_step2%in%sel_cls,]
#creaste psudeobulk groupings
pd_ATAC <- c()

for(x in sel_cls){
  del <- plot.data[plot.data$Cluster_step2==x,]
  cellsx <- row.names(del)
  nsplits <- floor(length(cellsx)/ 50)
  remainder <- length(cellsx)-(nsplits*50)
  sizes <- rep(50,nsplits)
  sizes[1] <- sizes[1] + remainder
  arr <- sample(cellsx) ## this will randomize it, else PSB are created by serial
  indices <- split(arr, rep(1:nsplits, times = sizes))
  del$Split <- "ND"
  for(i in 1:length(indices)){
    del[indices[[i]],"Split"] <- paste0("a",i)
  }
  del$PSB <- paste0(del$Cluster_step2,"_",del$Split)
  pd_ATAC <- rbind(pd_ATAC,del)
  rm(del,cellsx,nsplits,remainder,sizes,arr,indices)
}
rm(x)
table(row.names(plot.data) %in% row.names(pd_ATAC))
pd_ATAC <- pd_ATAC[row.names(plot.data),]

cells <- row.names(pd_ATAC)

##
peaks <- unique(regCREs$Region)

## process subsetted peak matrix------------
proj <- loadArchRProject(path = paste0(DIR_ATAC,"comATACx/"))

##extract PeakMatrix for selected cells
proj <- proj[cells,]
PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
rm(proj)

## get Matrix and subset to peaks
region_names <- as.data.frame(PeakMatrix@rowRanges)
row.names(region_names) <- paste0(region_names$seqnames,":",region_names$start,"-",region_names$end)
head(region_names)

PM <- assay(PeakMatrix)
row.names(PM) <- row.names(region_names)
PM[1:4,1:4]
PM <- PM[,cells] ###very important to make sure cell order is same

## AUC score for significance test-------------
ms <- AUCell_run(PM, GSEA, aucMaxRank=nrow(PM)*0.10, normAUC = FALSE,
                 BPPARAM=BiocParallel::MulticoreParam(9))
scr <- t(assay(ms))
save(plot.data,scr, file ="PNana/pdAUC_HOXCRE_GCPN.RData") #HOX GRNs
save(plot.data,scr, file ="PNana/pdAUC_TCFSMADCRE_GCPN.RData") #other GRNs
save(plot.data,scr, file ="PNana/pdAUC_NonHOXCRE4HOXtar_GCPN.RData") #CREs associated with HOX targets where other TFs binds

save(plot.data,scr, file ="PNana/pdAUC_DEGt100CRE_GCPN.RData") #top 100 gene cres, only focusing on CREs selected by SCENIC+

range_01 <- function(x){
  mx <- 1
  mn <- 0
  a <- (mx-mn)/(max(x)-min(x))
  b <- mx-(a*max(x))
  d <- round((a*x)+b, digits = 3)
  return(d)
  
}
## significance test from pairwise comparisons
marks <- pairwiseWilcox(t(scr), groups=plot.data$Cluster_step2,direction="up")
sets <- data.frame(A=marks$pairs$first,B=marks$pairs$second)
dd <- marks$statistics
names(dd) <- paste0(sets$A," vs ",sets$B)
dd

pval_df <- c()

for(x in sel_cls){
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

pd$Cluster <- factor(pd$Cluster, levels = sel_cls)

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

