## Subsetting the NMF factors obtained for entire hindbrain 
## to cells belonging to one of the 16 groups in lv1 clusters
## and increasing clustering resolution to obtain lv2 cluster
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")


#function for  return nearest neighbour less than radius and max neighbors
suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(Matrix)
})
#######Processing lineage clusters, getting NMF #############
##load data and process ----
sel_ann <- commandArgs(trailing=TRUE) ## one of the 16 lv1 groups
load("NMF_HB1708.RData")
load("plotdataHBv2.RData") ##fixed with old and CB annotation
keep <- plot.data$Annotation %in% sel_ann

plot.data <- plot.data[keep,]
rd_nmf <- rd_nmf[keep,]

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
knn = sk.NearestNeighbors(n_neighbors=6, radius=0.5,metric="precomputed", n_jobs=3) ##radius 0.5 
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
part = la.find_partition(g_rnn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=1000, seed=12)


k6 = part.modularity
#clusters
k7 = np.array(part.membership)

exit
k10 <- py$k7

print(py$k6)
#assigning clustering-
plot.data$RNN_cl <- paste0("progen_",k10)

save(plot.data, file = paste0("plotdata_",sel_ann,".RData"))
q()
