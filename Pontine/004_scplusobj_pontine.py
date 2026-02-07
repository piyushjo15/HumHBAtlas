
import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
import pickle

import numpy as np
import pybiomart as pbm

# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')


#Define dir
print("Start Analysis for Pontine Nuceli:")

projDir='/home/SCENICOut/PN/'
os.chdir(projDir)
tmpDir='/home/TMP/'


#Create Scenic object #############
# Load functions
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *

#load scRNA and scATAC data
adata = sc.read_h5ad(os.path.join(projDir,'adata_RNA.h5ad'))
cistopic_obj = dill.load(open(os.path.join(projDir, 'Cistopic/PN_cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(projDir, 'Cistopic/motifs/menr.pkl'),'rb'))

adata.raw = adata

#Create SCENIC+ object - for Non-Multiome data
scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr,
        multi_ome_mode = False,
        key_to_group_by = 'Cluster_step2',
        nr_cells_per_metacells = 5)


#scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj


print(f"The cell lines for which we have scRNA-seq data are:\t{', '.join(set(adata.obs['Cluster_step2']) - set(['-']))}")
print(f"The cell lines for which we have scATAC-seq data are:\t{', '.join(set(cistopic_obj.cell_data['Cluster_step2']))}")
print(f"The cell lines for which we have both:\t{', '.join(set(cistopic_obj.cell_data['Cluster_step2']) & set(adata.obs['Cluster_step2']))}")


# Save
if not os.path.exists(os.path.join(projDir, 'scenicplus')):
  os.makedirs(os.path.join(projDir, 'scenicplus'))


pickle.dump(scplus_obj, open(os.path.join(projDir, 'scenicplus/PN_scplus_obj.pkl'), 'wb'))

print("Saved scplus object")
