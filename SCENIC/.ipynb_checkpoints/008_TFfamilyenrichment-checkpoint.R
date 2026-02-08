
## This script perform analaysis for TF-family distribution

library(tidyverse)
library(ggplot2)
setwd("~/SCENIC")
load("AllGRNs.RData")

grns <- list()
cls <- unique(all_tfs$Class)
for(x in cls){
  grns[[x]] <- all_tfs[all_tfs$Class==x,"TF"]
}
lengths(grns)
library(UpSetR)


grns_input <- fromList(grns)

# make the upset plot
p <- upset(grns_input,
      nsets = 16,
      nintersects = NA,
      sets.x.label = "Elements",
      order.by = "freq")


tffam <- read.delim("tffamily.txt",row.names = 1) ##classfication of identified TFs into families
head(tffam)
head(all_tfs)
tfs <- unique(all_tfs$TF)
table(tfs %in% row.names(tffam))
all_tfs$Fam <- "ND"
for(x in tfs){
  all_tfs[all_tfs$TF==x,"Fam"] <- tffam[x,"Family"]
}


dx <- data.frame(table(tffam$Family))
dx <- dx[order(dx$Freq,decreasing = TRUE),]
dx$Var1 <- factor(dx$Var1, levels = dx$Var1)
colx <- colorRampPalette(pal_simpsons()(14))(length(unique(all_tfs$Fam)))
names(colx) <- unique(all_tfs$Fam)
p1 <- ggplot(dx,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat = "identity")+
scale_fill_manual(values=colx)+
  theme_bw()+
theme(axis.text.y=element_text(size=14,color="black"),
      axis.text.x=element_text(size=14,color="black",angle=45,hjust=1),
      legend.position="None")

fams <- data.frame(table(all_tfs$Class,all_tfs$Fam))
fams <- fams[fams$Freq!=0,]
fams$Var2 <- factor(fams$Var2,levels = dx$Var1)
cls <- c("Progenitors","Neuroblasts_I","Neuroblasts_II","Neuroblasts_III",
         "GCUBC","GC","CBINT","Purkinje","Neurons_I","Neurons_II",
         "Glioblasts","PreO","MOG","MixedGlia","MES_VAS","Immune")
all_tfs$Class <- factor(all_tfs$Class,levels = cls)

head(fams)
vx <- read.delim("Classcol.txt")
colx <- vx$Color
names(colx) <- vx$Class
library(ggsci)
colx <- colorRampPalette(pal_simpsons()(14))(length(unique(all_tfs$Fam)))
names(colx) <- unique(all_tfs$Fam)
p2<- ggplot(all_tfs,aes(x=Class, fill = Fam))+
  geom_bar(position = "fill")+
  scale_fill_manual(values=colx)+
  theme_bw()+
 theme(axis.line = element_line(color="black",linewidth=.8),
        axis.ticks = element_line(color="black",linewidth=.8),
        axis.text.x = element_text(color="black",angle=45,hjust=1,,size=14),
        axis.text.y=element_text(color="black",size=14),
        axis.title=element_blank(),
      legend.position="None")

gridExtra::grid.arrange(p1,p2)

# test for overrepresentation
## chi-square
tbl <- table(all_tfs$Class, all_tfs$Fam)
set.seed(123)
res <- chisq.test(tbl)
res$stdres[ abs(res$stdres) > 2 ] 
stdres <- res$stdres

sig <- which(abs(stdres) > 2, arr.ind = TRUE)
df <- data.frame(
  Class = rownames(stdres)[sig[,1]],
  Fam   = colnames(stdres)[sig[,2]],
  Resid = stdres[sig]
)

df_res <- as.data.frame(as.table(res$stdres))
colnames(df_res) <- c("Class", "Fam", "StdResid")
# 3. compute significance threshold (|resid| > 2 ≈ p < 0.05)
df_res$Significant <- abs(df_res$StdResid) > 2
# 4. plot
df_res$Fam <- factor(df_res$Fam, levels =dx$Var1)
df_res$Class <- factor(df_res$Class, levels =cls)
df_res$StdResid <- pmax(pmin(df_res$StdResid,2),-2)


# 4. plot

p <- ggplot(df_res, aes(y = Fam, x = Class)) +
  geom_point(aes(size = abs(StdResid), color = StdResid), alpha =.8) +
  geom_point(
    data = subset(df_res, Significant),
    aes(size = abs(StdResid)),
    shape = 21, stroke = .5, fill = NA, color = "black"
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0) +
  scale_size(range = c(1, 6)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color="black"),
    axis.text.y = element_text(size = 14, color="black"),
    axis.line = element_line(color="black", linewidth=.8),
    axis.ticks = element_line(color="black", linewidth=.8),
    axis.title = element_blank()
  )


##negative binomial
fams <- colnames(tbl)
df_counts <- as.data.frame(as.table(tbl))
df_counts <- df_counts[df_counts$Freq!=0,]
colnames(df_counts) <- c("Class","Family","Freq")


library(MASS)
m_add <- glm.nb(Freq ~ Class + Family, data = df_counts)

df_counts$stdres <- residuals(m_add, type = "pearson")
df_counts$stdres_capped <- pmax(pmin(df_counts$stdres, 2), -2)

# Flag significant cells (|z| > 2 ≈ p < 0.05)
df_counts$signif <- ifelse(abs(df_counts$stdres) > 2,
                           ifelse(df_counts$stdres > 0, "Over", "Under"),
                           "NS")
res_mat <- reshape2::acast(df_counts, Class ~ Family, value.var = "stdres_capped", fill = 0)
df_counts$Sel <- "No"
df_counts[df_counts$signif!="NS","Sel"] <- "Yes"

# Plot
df_counts$Family <- factor(df_counts$Family, levels = dx$Var1)
df_counts$Class <- factor(df_counts$Class, levels = cls)


p3 <- ggplot(df_counts, aes(Family, Class,color = stdres_capped,size  = abs(stdres_capped))) +
  geom_point() +
  geom_point(data = subset(df_counts, Sel == "Yes"),
             shape = 21, fill = NA, color = "black", size = 6, stroke = .5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        limits = c(-2, 2)) +
 theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color="black"),
    axis.text.y = element_text(size = 14, color="black"),
    axis.line = element_blank(),
    axis.ticks = element_line(color="black", linewidth=1),
    axis.title = element_blank(),
      panel.border=element_rect(color="black",linewidth=1))+
  coord_flip()
p3

