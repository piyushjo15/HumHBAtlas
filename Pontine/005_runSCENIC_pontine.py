import scenicplus
scenicplus.__version__

import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
import pickle

# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')


#Define dir
print("Start Analysis for Pontine Nuceli:")
projDir='/home/SCENICOut/PN/'
os.chdir(projDir)
tmpDir='/home/TMP/'

#load Scenicplus object
infile = open(os.path.join(projDir, 'scenicplus/PN_scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()


#Generate Cistromes
from scenicplus.cistromes import *

import time
start_time = time.time()
merge_cistromes(scplus_obj)
time = time.time()-start_time
print(time/60)

##### Infer Enhancer - Gene Relationships #######
#Get Search Space
biomart_host = "http://feb2021.archive.ensembl.org/" #paste correct biomart_host (defined in 004_create_scplus_obj.py script)

from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, RF_KWARGS
get_search_space(scplus_obj,
                 biomart_host = biomart_host,
                 species = 'hsapiens',
                 assembly = 'hg38',
                 upstream = [1000, 150000],
                 downstream = [1000, 150000])

#Enhancer to Gene Models using correlation = random forest or Gradient Boosting Machine

calculate_regions_to_genes_relationships(scplus_obj,
                                         ray_n_cpu = 2,
                                         _temp_dir = tmpDir,
                                         importance_scoring_method = 'RF',
                                         importance_scoring_kwargs = {'n_jobs': 3, 'max_features': 0.3, 'n_estimators': 1000})

# Save
pickle.dump(scplus_obj,
            open(os.path.join(projDir,'scenicplus/PN_scplus_obj.pkl'), 'wb'))

print("Saved scenicplus object after calculating regions-genes reslationships for pontine")


print("Add TF2G adjancecies to: ")
#load TSV file

file_path =os.path.join(projDir, 'corPN.tsv')
tf2g = pd.read_csv(file_path, sep='\t')
#remvoe NA from 'rho' column
tf2g_cleaned = tf2g.dropna(subset=['rho'])
# #save as TSV
output_file_path = os.path.join(projDir, 'corPN_cleaned.tsv')
tf2g_cleaned.to_csv(output_file_path, sep='\t', index=False)

from scenicplus.TF_to_gene import *
load_TF2G_adj_from_file(scplus_obj,
                        f_adj =output_file_path,
                        inplace = True,
                        key= 'TF2G_adj')

print('Checking if TF2G has been added to scplus_obj')
print(scplus_obj.uns.keys())

# Save
pickle.dump(scplus_obj,
            open(os.path.join(projDir,'scenicplus/PN_scplus_obj.pkl'), 'wb'))

print("Saved scenicplus object after adding TF2G adjacencies for ")

############# Build eGRNs ##############

# Load functions
from scenicplus.grn_builder.gsea_approach import build_grn


build_grn(scplus_obj,
         min_target_genes = 10,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.05,
         keep_extended_motif_annot = False, # False v2
         keep_only_activating=True, # True v2
         merge_eRegulons = True,
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance_x_rho',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = 4,
         _temp_dir = tmpDir)

print("Finished Build GRN")
# Save
#to access the eGRNs
with open(os.path.join(projDir, 'scenicplus/PN_scplus_obj.pkl'), 'wb') as f:
  dill.dump(scplus_obj, f)

print("saved after bulding eGRNs for Pontine:")

from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')
print("check unstructured data in scplus object after adding eRegulon metadata")
scplus_obj.uns.keys()

directory_path = os.path.join(projDir,'metadata/')

# Check if the directory already exists
if not os.path.exists(directory_path):
    os.makedirs(directory_path)


############# save metadata #####################

metadata = scplus_obj.uns['eRegulon_metadata']
metadata.to_csv(projDir + '/metadata/PN_eRegulon_metadata.csv',
                float_format='%.3f',index=False, sep=';')


print("Saved eRegulon metadata for Pontine:" )



