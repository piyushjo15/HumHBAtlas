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
## performic pairwise wilcoxon test for gene score activity
## AUC bubbple plot using module score-------------
### findinding modscore
load("pdPNGCUBC_GRNs_MS.RData")
#load("pdPNGC_DEGt100PNGC_MS.RData")
cl_lv <- c("Neuroblasts_A:A1_AES","Glioblasts:Progenitor_AES",
           "GCUBC:GCP_EM","GCUBC:GCP_PN")
plot.data <- pd[pd$Cluster%in%cl_lv,]
grns <- colnames(scr)
hox <- grep("HOX",grns,value = TRUE)

non_hox <- c("TCF7L1","TCF7L2","POU3F2","SMAD3","SMAD2")
non_hox <- paste0("PN:",non_hox)

scr <- scr[row.names(plot.data),hox] ##or non hox
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


