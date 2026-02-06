## Subsetting the NMF factors obtained for entire hindbrain 
## to cells belonging to one of the 16 groups in lv1 clusters
## and increasing clustering resolution to obtain lv2 cluster
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "~/.conda/envs/pjexv4/bin/python3.7")
use_python("~/.conda/envs/pjexv4/bin/python3.7")

#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(Matrix)
})
#######Processing lineage clusters, getting NMF #############
##load data and process ----
load("lgcountsHB.RData")

load("NMF_HB1708.RData")
load("plotdata_NNDCS_50k_180823fix.RData") ##fixed with old and CB annotation
keep <- plot.data$Annotation %in% "Progenitors"

plot.data <- plot.data[keep,]
rd_nmf <- rd_nmf[keep,]
mdt <- mdt[,row.names(plot.data)]
save(mdt, file="lgcounts_Progenitors.RData")
rm(mdt)
#q()
set.seed(787)
del <- uwot::umap(rd_nmf, n_neighbors=30, min_dist=0.3, metric="cosine")
plot.data$UMAP1 <- del[,1]
plot.data$UMAP2 <- del[,2]

############################# clustering###########################

repl_python()
import numpy as np
import igraph as ig
import leidenalg as la
import networkx as nx
import scipy.sparse as ss
import sklearn.neighbors as sk
import scipy.spatial.distance as ssd
import sklearn.metrics as skm
import random
# #using scipy cist for computing ditance matrix
X=r.rd_nmf
ax = ssd.cdist(X,X, metric="jensenshannon")
##using sklearn neighbors to genrate neighbors graph using precomputed sitance matrix
knn = sk.NearestNeighbors(n_neighbors=6, radius=0.5,metric="precomputed", n_jobs=3) ##radius 0.5 for JS

# # ##using sklearn neighbors to genrate neighbors graph using cosine metric
# knn = sk.NearestNeighbors(n_neighbors=6, radius=0.3,metric="cosine", n_jobs=3) ##radius 0.25 for cosine
# ax = r.rd_nmf
##fititng
knn.fit(ax)

##KNN
knnm = knn.kneighbors_graph(ax)
knnm= knnm.toarray()
#RNN
rnnm = knn.radius_neighbors_graph(ax)
rnnm= rnnm.toarray()

c1= np.multiply(rnnm,knnm) ##I have checked it and this works
##type(c1) ##for class

del(ax,knn,knnm,rnnm)
## creating a networkx graph
nxg = nx.from_numpy_array(c1)

#converting to igraph
g_rnn = ig.Graph(len(nxg), list(zip(*list(zip(*nx.to_edgelist(nxg)))[:2])), directed=False)
del(nxg)

#getting edgelist
ela = g_rnn.get_edgelist()
## gettign weight list
wj = g_rnn.similarity_jaccard(pairs=ela) #for edge list pairs
#getting clustering
#part = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=1000, seed=12)
#same as above but with resolution parameter
#part = la.find_partition(g_rnn, la.RBConfigurationVertexPartition, resolution_parameter = 3, weights=wj, n_iterations= -1, max_comm_size=1000, seed=12)

#another approach
opt = la.Optimiser()
opt.consider_empty_community = True
opt.max_comm_size = 1000
opt.set_rng_seed(123)
part = la.RBConfigurationVertexPartition(g_rnn,  resolution_parameter = 2, weights=wj)
opt.optimise_partition(part, n_iterations=-1)

k6 = part.modularity
#clusters
k7 = np.array(part.membership)

exit
k10 <- py$k7

print(py$k6)
#assigning clustering-
plot.data$RNN_cl <- paste0("progen_",k10)

save(plot.data, file = "plotdata_Progenitors_RNN060124.RData" )
q()
###----------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggsci)
  library(ggrepel)
  library(gridExtra)
})
load("plotdata_Progenitors.RData")
d <- length(unique(plot.data$Batch))
#plot.data$NND_cl <- factor(plot.data$NND_cl, levels = paste0("NNDJS_",c(0:(d-1))))
label.d = plot.data %>% group_by(Batch) %>% 
  select(UMAP1, UMAP2) %>% summarize_all(median)
#for plot
ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=Batch))+
  geom_point(size=0.5)+
  scale_color_manual(values =colorRampPalette(pal_igv()(51))(d))+
  theme_classic()+
  geom_label_repel(aes(label = Batch),size = 2.5, data = label.d, show.legend = FALSE)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank(),
        legend.position = "none")



st <- read.delim("Stagecol.txt", row.names = 1)
plot.data$Stage<-factor(plot.data$Stage,levels=row.names(st))

valuesx <-  st$Color
ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=Stage))+
  geom_point(size=0.5)+
  scale_color_manual(values=valuesx)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank())
ggplot(plot.data, aes(x=UMAP1, y=UMAP2, color=Region))+
  geom_point(size=0.5)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title=element_blank())
