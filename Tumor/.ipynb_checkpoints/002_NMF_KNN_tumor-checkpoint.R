##Performing NMF factorization followed by KNN and clustering for PA/DIPG/PFA samples
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")


suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(Matrix)
  library(foreach)
})
DIR_OUT="~/DIPG/"
#DIR_OUT="~/PFA/"
#DIR_OUT="~/PA/"

setwd(DIR_OUT)
#
# ## 1. load multibatchnormed cosnormed  merged----
# load("cosnormDIPG.RData")
# mdt <- t(com.sce_cnnc)
# rm(com.sce_cnnc)
# load("plotdataDIPG10.RData")
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
# save(H, rd_nmf, file = "NMF_DIPG.RData" )
# #q()
# ## run UMAP
# set.seed(787)
# del <- uwot::umap(rd_nmf, n_neighbors=30, min_dist=0.3)
# 
# plot.data$UMAP1 <- del[,1]
# plot.data$UMAP2 <- del[,2]
# 
# save(plot.data, file = "plotdataDIPG10v2.RData")
# #rm(rd_nmf, H)
# 
## python code ----
load("NMF_DIPG.RData")
load("plotdataDIPG10v2.RData")
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
save(plot.data, file="plotdataDIPG_KNN.RData")
print(py$k6)
q()
##post plotdataDIPG_KNN.RData, KNN_cl were annotated into OPC, Astro, Microglia, and other identities into plotdataDIPG_ANN.RData, and same for PA and GJPFA