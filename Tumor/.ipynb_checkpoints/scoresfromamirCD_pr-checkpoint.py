import os
import pickle
import pandas as pd
import networkx as nx
import sys
import numpy as np
## Per sample
SAM=sys.argv[1]
proj_DIR="/home/p541i/MBsnANA/Diemana/scepro/RECORDR/"
os.chdir(proj_DIR)
print("Performing Word2Vec analysis for sample "+SAM)

partx='a'
## load the cosine drift scores
#with open(f"PA_cosines_pr_{partx}.pkl", "rb") as f:
#  df=pickle.load(f)
with open(f"DMG_cosines_pr_{partx}.pkl", "rb") as f:
  df=pickle.load(f)

#with open(f"PFA_cosines_pr_{partx}.pkl", "rb") as f:
#  df=pickle.load(f)

df=df.T ##trasnpose
## this is cosine similairty, 1-cosine distance
cosine_dict = df[SAM].to_dict() ##obtain dict of score for a sample


#load graph


with open(f"Normalv2/gene_network_nor_pr_{partx}.pkl", "rb") as f:
  G_nor=pickle.load(f)
#with open(f"Normal_PFA/gene_network_nor_pr_{partx}.pkl", "rb") as f:
#  G_nor=pickle.load(f)
with open(f"{SAM}/gene_network_{SAM}_pr.pkl", "rb") as f:
    G_tum=pickle.load(f)

# Set of all genes (nodes) in relapse graph
genes_tum = set(G_tum.nodes())
genes_nor = set(G_nor.nodes())


## add misisng tumor gene cosine drift score
for gene in genes_tum:
    if gene not in cosine_dict:
        cosine_dict[gene] = 1.2

# Create a cosine lookup dictionary for index gene scores

# Initialize results
results = []
from concurrent.futures import ProcessPoolExecutor
all_neighbors = {node: set(G_tum.neighbors(node)) for node in G_tum.nodes()}

def score_gene(gene):
    neighbors_tum = all_neighbors.get(gene, set())
    neighbors_nor = set(G_nor.neighbors(gene)) if gene in G_nor else set()
    unique_neighbors = neighbors_tum - neighbors_nor
    neighbor_score = len(unique_neighbors) / len(genes_tum)

    reachable_genes = set(neighbors_tum)
    second_hop = set.union(*(all_neighbors.get(n, set()) for n in neighbors_tum))
    reachable_genes.update(second_hop)
    graph_reach_score = len(reachable_genes) / len(genes_tum)

    valid_neighbors = [n for n in neighbors_tum if n in cosine_dict and -1 <= cosine_dict[n] <= 1]
    if valid_neighbors:
        mean_cosine = np.mean([cosine_dict[n] for n in valid_neighbors]) ## cosine similairty
        context_drift_score = 1 - mean_cosine ## here I need cosine distance, because I want to reward differences and not similairty
    else:
        context_drift_score = 0

    index_gene_score = 1- cosine_dict[gene] if gene in cosine_dict and -1 <= cosine_dict[gene] <= 1 else 1 ## here I need cosine distance

    total_score = neighbor_score + graph_reach_score + context_drift_score + index_gene_score

    return {
        'gene': gene,
        'neighbor_score': neighbor_score,
        'graph_reach_score': graph_reach_score,
        'neighborhood_context_drift': context_drift_score,
        'index_gene_context_drift': index_gene_score,
        'total_score': total_score
    }

with ProcessPoolExecutor() as executor:
    results = list(executor.map(score_gene, genes_tum))


# Create final DataFrame
df_scores = pd.DataFrame(results)
df_scores.sort_values('total_score', ascending=False, inplace=True)

df_scores.to_csv(f"{SAM}/gene_context_drift_scores_nor_{SAM}_pr_amir_{partx}.csv", index=False)
