import os
import pickle
import numpy as np
import pandas

dataDir='/home/DeepHB/'
os.chdir(dataDir)

def parse_jaspar_meme(filename):
    motifs = {}
    with open(filename, 'r') as f:
        lines = f.readlines()

    motif_id = None
    matrix = []
    reading_matrix = False

    for line in lines:
        line = line.strip()
        
        if line.startswith('MOTIF'):
            # Save previous motif if it exists
            if motif_id and matrix:
                motifs[motif_id] = matrix
            # Start new motif
            motif_id = line.split()[1]
            matrix = []
            reading_matrix = False
        
        elif line.startswith('letter-probability matrix:'):
            reading_matrix = True
            continue

        elif reading_matrix:
            if line == '' or line.startswith('MOTIF') or line.startswith('URL'):
                reading_matrix = False
                if motif_id and matrix:
                    motifs[motif_id] = matrix
                matrix = []
            else:
                values = list(map(float, line.strip().split()))
                matrix.append(values)

    # Don't forget the last motif
    if motif_id and matrix:
        motifs[motif_id] = matrix

    return motifs


motifs_dict = parse_jaspar_meme("extrafiles/JASPAR2024_root.meme")
np.array(motifs_dict["cluster_001"])

# Save to pickle
with open("extrafiles/JASPAR2024_cluster_motifs_dict.pkl", "wb") as f:
    pickle.dump(motifs_dict, f)
