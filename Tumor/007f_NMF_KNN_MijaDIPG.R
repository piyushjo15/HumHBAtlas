##Performing NMF factorization followed by KNN and clustering for PA samples
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/p541i/MBsnANA/conda_env/pyforR/bin/python")
use_python("/home/p541i/MBsnANA/conda_env/pyforR/bin/python", required = TRUE)

suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(Matrix)
  library(foreach)
})
DIR_OUT="~/MBsnANA/HGGana/Mija_DIPG/"
setwd(DIR_OUT)
#
# ## 1. load multibatchnormed cosnormed  merged----
# load("cosnormMijaDIPG10x.RData")
# mdt <- t(com.sce_cnnc)
# rm(com.sce_cnnc)
# load("plotdataMijaDIPG10x170625.RData")
# 
# ## 2. running NMF using sklearn----
# repl_python()
# import sklearn.decomposition as sk
# import numpy
# import random
# # set seed and run NMF
# # KL solver not working.
# model = sk.NMF(n_components=25, init='nndsvda', max_iter=10000)
# 
# random.seed(456)
# w = model.fit_transform(r.mdt)
# h = model.components_
# exit
# 
# rd_nmf= py$w
# H = py$h
# #row.names(rd_nmf) <- row.names(plot.data)
# #colnames(rd_nmf) <- row.names(H) <- paste0("NMF_",seq(dim(rd_nmf)[2]))
# #row.names(H) <- colnames(mdt)
# save(H, rd_nmf, file = "NMF_MijaDIPG.RData" )
# #q()
# ## run UMAP
# set.seed(787)
# del <- uwot::umap(rd_nmf, n_neighbors=30, min_dist=0.3)
# 
# plot.data$UMAP1 <- del[,1]
# plot.data$UMAP2 <- del[,2]
# 
# save(plot.data, file = "plotdataMijaDIPG10x170625v2.RData")
# #rm(rd_nmf, H)
# 
## python code ----
load("NMF_MijaDIPG.RData")
load("plotdataMijaDIPG10x170625v2.RData")
repl_python()
import numpy as np
import umap
import random
import igraph as ig
import leidenalg as la
import networkx as nx
import scipy.sparse as ss
import sklearn.neighbors as sk

wn=r.rd_nmf

##using NMF factors to identify nearest neighbors using sklearn neighbors graph
## This uses data to give a neighborhood graph
random.seed(455)
a = sk.kneighbors_graph(wn,n_neighbors=26, metric="cosine", n_jobs=3, include_self=True)

## creating a networkx graph
nxg = nx.from_scipy_sparse_array(a) ##convert to networkx grapg
#converting to igraph
g_rnn = ig.Graph(len(nxg), list(zip(*list(zip(*nx.to_edgelist(nxg)))[:2])), directed=False)
del(nxg,a)

#getting edgelist
ela = g_rnn.get_edgelist()
## gettign weight list
wj = g_rnn.similarity_jaccard(pairs=ela) #for edge list pairs
#getting clustering
part = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=2000, seed=12)

#another approach for resolution ## not doing this for PA
# opt = la.Optimiser()
# opt.consider_empty_community = True
# opt.max_comm_size = 1000
# opt.set_rng_seed(123)
# part = la.RBConfigurationVertexPartition(g_rnn,  resolution_parameter = 0.3, weights=wj)
# opt.optimise_partition(part, n_iterations=-1)

k6 = part.modularity
#clusters
k7 = np.array(part.membership)

exit
k10 <- py$k7

print(py$k6)
#assigning clustering-
plot.data$KNN_cl <- paste0("KNN_",k10)
save(plot.data, file="plotdataMijaDIPG10x_KNN_1706.RData")
print(py$k6)
q()
## check -----
suppressPackageStartupMessages({
  library(ggrepel)
  library(ggsci)
  library(dplyr)
  library(tidyr)
})
load("plotdataMijaDIPG10x_ANN.RData")

ann <- read.delim("Mija_KNN.txt",row.names = 1)
cl_lv <- row.names(ann)
head(ann)
anns <- row.names(ann)
plot.data$ANN1 <- plot.data$ANN2 <- "ND"
for(x in anns){
  plot.data[plot.data$KNN_cl==x,"ANN1"] <- ann[x,"ANN1"]
  plot.data[plot.data$KNN_cl==x,"ANN2"] <- ann[x,"ANN2"]
  plot.data[plot.data$KNN_cl==x,"ANN3"] <- ann[x,"ANN3"]
}
table(plot.data$ANN1)
table(plot.data$ANN2)
table(plot.data$ANN3,plot.data$ANN1)

d <- length(unique(plot.data$Batch))
label.d = plot.data %>% group_by(Batch) %>%
  select(UMAP1, UMAP2) %>% summarize_all(median)

# #for plot
p <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=Batch))+
  geom_point(size=1)+
  #scale_color_manual(values =valuesx)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Batch),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")

pdf("Figures/UMAP_Batch.pdf", width=8, height=8,pointsize = 10)
print(p)
dev.off()

d <- length(unique(plot.data$KNN_cl))
label.d = plot.data %>% group_by(KNN_cl) %>%
  select(UMAP1, UMAP2) %>% summarize_all(median)

# #for plot
p <- ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=KNN_cl))+
  geom_point(size=1)+
  scale_color_manual(values =colorRampPalette(pal_simpsons()(16))(d))+
  theme_classic()+
  geom_label_repel(aes(label = KNN_cl),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "None")
tiff("Figures/UMAP_Batch.tiff", units="in", width=8, height=8, res=150)
print(p)
dev.off()
save(plot.data, file="plotdataPAnew10x_KNN_1008.RData")


values <- read.delim("~/MBsnANA/HGGana/Scripts/RNA/glial_ann_col.txt",row.names = 1)
cls <- row.names(values)
plot.data$ANN <- plot.data$SVM_pri_old
k <- grep("Immune",plot.data$ANN)
plot.data[k,"ANN"] <- "Immune"
keep <- plot.data$ANN %in% cls
table(keep)
plot.data[!keep,"ANN"] <- "Others"
d <- sort(unique(plot.data$ANN))
valuesx <- values[d,1]
table(is.na(valuesx))

p <- ggplot(plot.data, aes(x=KNN_cl, fill=ANN))+
  scale_fill_manual(values =valuesx)+
  geom_bar(position = "fill")+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth=1),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.x = element_text(face="bold",colour = "black",size=10, angle = 90,vjust=0.5, hjust=0),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title=element_blank())+
  scale_y_continuous(breaks =c(0, 0.5, 1))+
  facet_wrap(~Batch,ncol=2, scales = "free_x")

tiff("JSH3HGG_ANN2.tiff", units="in", width=8, height=4, res=300)
print(p)
dev.off()

ggplot(plot.data, aes(x=KNN_cl, y=decontX_contamination))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth=1),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.x = element_text(face="bold",colour = "black",size=10, angle = 90,vjust=0.5, hjust=0),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title=element_blank())
d <- length(unique(plot.data$Batch))

ggplot(plot.data, aes(x=KNN_cl, fill=Batch))+
  geom_bar(position = "fill")+
  scale_fill_manual(values =c("red","blue","green","yellow","grey"))+
  theme_classic()+
  theme(axis.line = element_line(colour = 'black',linewidth=1),
        axis.ticks = element_line(colour = 'black',linewidth=1),
        axis.text.x = element_text(face="bold",colour = "black",size=10, angle = 90,vjust=0.5, hjust=0),
        axis.text.y = element_text(face="bold",colour = "black",size=8),
        axis.title=element_blank())
