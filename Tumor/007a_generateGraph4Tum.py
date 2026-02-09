import os
import sys
import networkx as nx
import pandas as pd
import numpy as np
import pickle
import community
from scipy.sparse import load_npz
from sklearn.metrics.pairwise import cosine_similarity
import leidenalg as la
import igraph as ig
os.environ['PYTHONHASHSEED'] = '42'
##generating a gene-gene community graph for each tumor, PA, DIPG or PFA
## Per sample
SAM=sys.argv[1]
proj_DIR="/home/RECORDR/"+SAM
os.chdir(proj_DIR)
print("Performing Node2Vec analysis for sample "+SAM)

## 1.1 Load tumor gene expression data
mat_tum = load_npz(f"genemat_"+SAM+"_pc.npz")
genes_tum = np.load(f"genes_pc.npy", allow_pickle=True)
print(mat_tum.shape)
print(genes_tum[0:5])
print(len(genes_tum))

## 1.2 Create correlation matrix and obtain node-edge list
arr = mat_tum.toarray()
rowvar = (arr.shape[0] == len(genes_tum))  # assume rows are genes iff rows==#genes
sim = np.corrcoef(arr, rowvar=rowvar)

#Clean & prep (make sure it's symmetric, drop self-edges)
sim = np.nan_to_num(sim, nan=0.0)
sim = 0.5 * (sim + sim.T)  # guard against tiny asymmetries
np.fill_diagonal(sim, 0.0)

# Build edge list from the upper triangle only with threshold
thr = 0.1
iu = np.triu_indices_from(sim, k=1)
vals = sim[iu]
mask = vals >= thr

df_tum = pd.DataFrame({
    "Node":   genes_tum[iu[0][mask]],
    "Target": genes_tum[iu[1][mask]],
    "correlation": vals[mask].astype(np.float32)
})

#  Save
with open(f"gene_edge_list_{SAM}.pkl", "wb") as f:
    pickle.dump(df_tum, f)

## 1.3 Create graph from adjacency list
## reading edge list
#df_tum= pd.read_pickle(f"gene_edge_list_{SAM}.pkl")

##created weighted Graph from adjacency list
# 1) Build weighted graph
G_tum = nx.from_pandas_edgelist(
    df_tum, 'Node', 'Target', edge_attr=['correlation'], create_using=nx.Graph()
)

for u, v, d in G_tum.edges(data=True):
    r = float(d['correlation'])
    # Soft-threshold:
    d['weight'] = r ** 6
    

## using community parition
part_tum = community.best_partition(G_tum, weight='weight',resolution=1,random_state=42)                      

nx.set_node_attributes(G_tum, part_tum, 'module')
# # Save filtered outputs
nx.write_gexf(G_tum, f"gene_network_{SAM}.gexf")
nx.write_edgelist(G_tum, f"gene_network_{SAM}.edg", delimiter='\t', data=['weight'])
with open(f"gene_network_{SAM}.pkl", "wb") as f:
    pickle.dump(G_tum, f)


# Save filtered communtiies
df_partition = pd.DataFrame(list(part_tum.items()), columns=["gene", "module"])
df_partition.to_csv(f"partition_{SAM}.csv", index=False)
# 
np.save(f"genes_fil_{SAM}.npy", np.array(df_partition.gene))

print("Finsihed generating fraph for sample "+SAM+"..")

exit()
