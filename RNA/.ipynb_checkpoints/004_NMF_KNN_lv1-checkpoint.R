##Performing NMF factorization of cosinenormed gene-expression data
## followed by RNN and clustering
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/home/loc/to/python/bin/python3.7")
use_python("/home/loc/to/python/bin/python3.7")

suppressPackageStartupMessages({
  library(batchelor)
  library(scran)
  library(scater)
  library(Matrix)
})

## load multibatchnormed cosnormed  merged----
load("cosnormHB.RData")
load("plotdataHB.RData")

dim(com.sce_cnnc)
mdt <- t(com.sce_cnnc) 
genes <- row.names(com.sce_cnnc)
rm(com.sce_cnnc)
repl_python()
import sklearn.decomposition as sk
import numpy
import random
# set seed and run NMF
# KL solver not working.
model = sk.NMF(n_components=100, init='nndsvda', max_iter=10000, random_state=0) 
random.seed(456)
w = model.fit_transform(r.mdt)
h = model.components_
exit

rd_nmf= py$w
H = py$h
save(H, rd_nmf, genes, file = "NMF_HB.RData" )
#q()
## run UMAP
set.seed(787)
del <- uwot::umap(rd_nmf, n_neighbors=50, min_dist=0.3, metric = "cosine")


plot.data$nfUMAP1 <- del[,1]
plot.data$nfUMAP2 <- del[,2]
#table(row.names(plot.data)==colnames(com.sce_rnnc))

save(plot.data, file = "plotdataHBv2.RData")

rm(rd_nmf, H)
#q()
###clustering -------
#load("NMF_HB.RData")
#load("plotdataHBv2.RData")
## python code ----
repl_python()
import pynndescent as nd
import numba
import numpy as np
import igraph as ig
import leidenalg as la
import networkx as nx
import scipy.sparse as ss
import random
import sklearn.neighbors as skn
#import time
##defining jensen shanon distance
#have verified with scipy distance

@numba.vectorize(fastmath=True)
def relative_entropy(x, y):
  if np.isnan(x) or np.isnan(y):
    return np.nan
  elif x > 0 and y > 0:
    return x * np.log(x / y)
  elif x == 0 and y >= 0:
    return 0.0
  else:
    return np.inf

@numba.njit(nogil=True, fastmath=True)
def jensen_shannon(p, q, base=None):
  p = p / np.sum(p, axis=0)
  q = q / np.sum(q, axis=0)
  m = (p + q) / 2.0
  left = relative_entropy(p, m)
  right = relative_entropy(q, m)
  js = np.sum(left, axis=0) + np.sum(right, axis=0)
  if base is not None:
    js /= np.log(base)
  return np.sqrt(js / 2.0)
  
#calling function to compile it, so it runs faster later
# w = r.rd_nmf
# x = w[0,0:10]
# y = w[1,0:10]
# jensen_shannon(x,y)
##oldand time consuming method
# ## Calling approx nearest neighbor function with jensen metric 
# random.seed(123)
# Y = nd.NNDescent(w, n_neighbors=50, metric=jensen_shannon)
# # getting indices of nighbors and node
# indices = Y.neighbor_graph[0]
# #generating a 0,1 neighbor graph
# x = w.shape[0]
# a = ss.csr_matrix((x,x), dtype=np.float32)
# for i in range(0,x):
#   a[i,indices[i,]] = 1

##newer way of running the matrix
random.seed(123)
Y = nd.NNDescent(w, n_neighbors=50, metric=jensen_shannon, random_state=101, low_memory=False, n_jobs=-1)
# getting indices of nighbors and node
indices = Y.neighbor_graph[0]
#generating a 0,1 neighbor graph
x = w.shape[0]
a = ss.lil_matrix((x,x), dtype=np.float32)

#start = time.time()
for i in np.arange(0,x):
  a[i,indices[i,]] = 1

#print(f'For loop: {time.time() - start} seconds')

## creating a networkx graph
nxg = nx.from_scipy_sparse_matrix(a)

#converting to igraph
g_nn = ig.Graph(len(nxg), list(zip(*list(zip(*nx.to_edgelist(nxg)))[:2])), directed=False)
del(nxg,a)
#getting edgelist
ela = g_nn.get_edgelist()
## gettign weight list
#wj <- g3$similarity_jaccard() #for all pairs
wj = g_nn.similarity_jaccard(pairs=ela) #for edge list pairs
#getting clustering
k5 = la.find_partition(g_nn, la.ModularityVertexPartition, weights=wj, n_iterations= -1, max_comm_size=4000, seed=12)

k6 = k5.modularity
#clusters
k7 = np.array(k5.membership)
exit

#back to r ----
#clusters 
cls <- py$k7

plot.data$NND_cl <- paste0("NNDJS_",cls)
save(plot.data, file="plotdataHBv2.RData")
print(py$k6)
q()
