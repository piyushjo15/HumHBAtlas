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

## I have already added Cistromes, performed region to gene linkages and added TF2adj in previous analysis
## just running buld eGRN. I am also not saving the changes in scenicplus object
## In the trial I played with top_n_regionTogenes_per_gene and top_n_regionTogenes_per_region
## making them (5,10) each.
## I also changes the NES_thr=1 which reduced the number of GRNs
## I also made chnages to keep_only_activating and keep_extended_motif_annot values
## the most important change is order_TFs_to_genes_by = 'importance_x_rho'

#Define dir
Cls = sys.argv[1]
Splt=sys.argv[2] ## 'A' or 'B'
print("Creating scenicplus object for Class: " + Cls+" and split: "+Splt)

DIROUT='/home/p541i/MBsnANA/HBana/SCENICana/SCENICOut/Class/'
adjDir=os.path.join(DIROUT,Cls)
projDir = os.path.join(DIROUT,Cls,'SCENIC',Cls+'_'+Splt)

tmpDir='/home/p541i/MBsnANA/TMP/'

#load Scenicplus object
infile = open(os.path.join(projDir,Cls+'_'+Splt+'_scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()


############# Build eGRNs ##############

# Load functions
from scenicplus.grn_builder.gsea_approach import build_grn


build_grn(scplus_obj,
         min_target_genes = 20,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.03,
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

from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')

############# save metadata #####################

metadata = scplus_obj.uns['eRegulon_metadata']
metadata.to_csv(projDir + '/GRNs/' +Cls+'_'+Splt+'_eRegulon_metadatav2.csv',
                float_format='%.3f',index=False, sep=';')


print("Saved eRegulon metadata for sample:" + Cls+'_'+Splt+' for version v2')
exit()


