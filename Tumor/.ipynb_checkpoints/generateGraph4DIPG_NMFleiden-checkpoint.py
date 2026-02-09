## In this version I am usign NMF factorization to obtain gene-to-gene connections
## using this to create the graph and gene communities as usual using leiden
## I am also making a weighted graph this time, which will be the cosine similarlirty 
## , hoping the weight impact the random walks.

import os
import sys
import networkx as nx
import pandas as pd
import numpy as np
import pickle
from scipy.sparse import load_npz
import leidenalg as la
import igraph as ig
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize
from scipy.spatial.distance import cosine

SAM=sys.argv[1]
proj_DIR="/home/p541i/MBsnANA/Diemana/scepro/RECORDR/"+SAM
os.chdir(proj_DIR)
print("Performing Node2Vec analysis for sample "+SAM)
## 1.1 Load gene expression data

##tumor
mat_tum = load_npz("genemat_"+SAM+"_pc.npz")
genes_tum = np.load("genes_pc.npy", allow_pickle=True)
print(mat_tum.shape)
print(genes_tum[0:5])
print(len(genes_tum))

## 1.2 Perform NMF, obtain H, normalize, calculate cosine sim, create edge list

nmf_mod = NMF(n_components=100, init='random', random_state=42, max_iter=1000)
W = nmf_mod.fit_transform(mat_tum.T)  # cells x factors
H = nmf_mod.components_  # factors x genes

np.save(f"H_{SAM}_r100.npy", H)

#H = np.load(f"H_{SAM}_r100.npy", allow_pickle=False)


from collections import defaultdict
from itertools import combinations

## connecting genes based on shared factor and correlation
top_n = 500
factor_gene_sets = []
for factor_row in H:
    top_indices = np.argsort(factor_row)[-top_n:]
    factor_gene_sets.append(set(top_indices))

co_occurrence = defaultdict(set)
for gene_set in factor_gene_sets:
    for i, j in combinations(gene_set, 2):
        co_occurrence[i].add(j)
        co_occurrence[j].add(i)

corr_mat = np.corrcoef(mat_tum.toarray())
#corr_mat = np.corrcoef(H.T)

# 3. Build edges with Pearson correlation ≥ 0.03
edges = []
corr_thres = 0.03

for i in co_occurrence:
    for j in co_occurrence[i]:
        if i < j:  # avoid duplicates
            corr = corr_mat[i, j]
            if corr >= corr_thres:
                edges.append((genes_tum[i], genes_tum[j], corr))


with open(f"gene_edge_list_{SAM}_nmf.pkl", "wb") as f:
    pickle.dump(edges, f)
## 1.3 Create graph from edge list

G_tum = nx.Graph()
G_tum.add_weighted_edges_from(edges)


#### leiden clustering to break up large communities
# Map node names to integers
node_names = list(G_tum.nodes())
name_to_id = {name: i for i, name in enumerate(node_names)}

# Create edges with integer IDs
edges_weighted = [(name_to_id[u], name_to_id[v], G_tum[u][v]['weight']) for u, v in G_tum.edges()]

# Build igraph graph
edges = [(u, v) for u, v, w in edges_weighted]
weights = [w for u, v, w in edges_weighted]

# Build graph with edges
g_nn = ig.Graph(edges=edges, directed=False)
# Add weights as edge attribute
g_nn.es['weight'] = weights
# Add original names as vertex attribute
g_nn.vs['name'] = node_names
k6 = la.find_partition(g_nn, la.ModularityVertexPartition,n_iterations =-1,
                      max_comm_size=1000,seed=42)
                      
## filter nodes
part_tum = {g_nn.vs[i]['name']: comm for i, comm in enumerate(k6.membership)}
# Save unfiltered
df_partition = pd.DataFrame(list(part_tum.items()), columns=["gene", "module"])
df_partition.to_csv(f"partition_{SAM}_nmf.csv", index=False)
del(df_partition)

nx.set_node_attributes(G_tum, part_tum, 'module')

from collections import Counter
# Filter small communities (size < 20)
min_size = 20
comm_sizes = Counter(part_tum.values())
nodes_to_keep = [node for node, comm in part_tum.items() if comm_sizes[comm] >= min_size]


# Filtered graph
G_filtered = G_tum.subgraph(nodes_to_keep).copy()
part_f = {node: part_tum[node] for node in nodes_to_keep}

# Save filtered outputs
nx.write_gexf(G_filtered, f"gene_network_{SAM}_nmf_fil.gexf")
nx.write_edgelist(G_filtered, f"gene_network_{SAM}_nmf_fil.edg", delimiter='\t', data=['weight'])
with open(f"gene_network_{SAM}_nmf_fil.pkl", "wb") as f:
    pickle.dump(G_filtered, f)

np.save("genes_nmf_fil.npy", np.array(nodes_to_keep))


# Save filtered communtiies
df_partition = pd.DataFrame(list(part_f.items()), columns=["gene", "module"])
df_partition.to_csv(f"partition_{SAM}_nmf_fil.csv", index=False)

del(df_partition)   
    

print("Finsihed generating fraph for sample "+SAM+"..")

exit()
