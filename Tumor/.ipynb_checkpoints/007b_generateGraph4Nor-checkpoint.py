## this sript generate GN for eubsetted glial lineage
## I have to grph as edge list for pecanpy
import os
import sys
import networkx as nx
import pandas as pd
import numpy as np
import pickle
import community
from scipy.sparse import load_npz
import leidenalg as la
import igraph as ig
from scipy.spatial.distance import cosine
os.environ['PYTHONHASHSEED'] = '42'

proj_DIR="/RECORDR/Normal/"

os.chdir(proj_DIR)
print("Performing Node2Vec analysis for Normal samples ")
partx='c'
## 1.1 Load gene expression data
mat_nor = load_npz(f"genemat_nor_pc_{partx}.npz")
genes_nor = np.load(f"genes_pc_{partx}.npy", allow_pickle=True)

print(mat_nor.shape)
print(genes_nor[0:5])
print(len(genes_nor))

## 1.2 Create correlation matrix and obtain node-edge list
arr = mat_nor.toarray()
rowvar = (arr.shape[0] == len(genes_nor))  # assume rows are genes iff rows==#genes
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

df_nor = pd.DataFrame({
    "Node":   genes_nor[iu[0][mask]],
    "Target": genes_nor[iu[1][mask]],
    "correlation": vals[mask].astype(np.float32)
})


# Save DataFrame as pickle
with open(f"gene_edge_list_nor_{partx}.pkl", "wb") as f:
  pickle.dump(df_nor, f)

#df_nor= pd.read_pickle(f"gene_edge_list_nor_{partx}.pkl")

##created weighted Graph from adjacency list
# 1) Build weighted graph
G_nor = nx.from_pandas_edgelist(
    df_nor, 'Node', 'Target', edge_attr=['correlation'], create_using=nx.Graph()
)

for u, v, d in G_nor.edges(data=True):
    r = float(d['correlation'])
    # Soft-threshold:
    d['weight'] = r ** 6
    

## using community parition
part_nor = community.best_partition(G_nor, weight='weight',resolution=1,random_state=42)                      
nx.set_node_attributes(G_nor, part_nor, 'module')

# Save filtered outputs
nx.write_gexf(G_nor, f"gene_network_nor_{partx}.gexf")
nx.write_edgelist(G_nor, f"gene_network_nor_{partx}.edg", delimiter='\t', data=['weight'])
with open(f"gene_network_nor_{partx}.pkl", "wb") as f:
    pickle.dump(G_nor, f)

# Save filtered communtiies
df_partition = pd.DataFrame(list(part_nor.items()), columns=["gene", "module"])
df_partition.to_csv(f"partition_nor_{partx}.csv", index=False)

# 
np.save(f"genes_fil_{partx}.npy", np.array(df_partition.gene))


print("Finsihed generating graph for Normal..")

exit()
