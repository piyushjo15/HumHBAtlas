
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

import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

## create SCENIC+ plus object (non-multiomic) for each class and split
#Define dir
Cls = sys.argv[1]
Splt=sys.argv[2] ## 'A' or 'B'
print("Creating scenicplus object for Class: " + Cls+" and split: "+Splt)

DIROUT='/home/SCENICOut/Class/'
projDir = os.path.join(DIROUT,Cls,'SCENIC',Cls+'_'+Splt)

tmpDir='/home/p541i/MBsnANA/TMP/'

#Create Scenic object #############
# Load functions
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *

#load scRNA and scATAC data
adata = sc.read_h5ad(os.path.join(DIROUT,Cls+'_'+Splt+'_RNA.h5ad'))
cistopic_obj = dill.load(open(os.path.join(projDir, Cls+'_'+Splt+'_cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(projDir, 'menr.pkl'),'rb'))

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


pickle.dump(scplus_obj, open(os.path.join(projDir,Cls+'_'+Splt+'_scplus_obj.pkl'), 'wb'))

print("Saved scplus object")

#Check with which Biomart host the gene names used in analysis match best
ensembl_version_dict = {'105': 'http://www.ensembl.org',
  '104': 'http://may2021.archive.ensembl.org/',
  '103': 'http://feb2021.archive.ensembl.org/',
  '102': 'http://nov2020.archive.ensembl.org/',
  '101': 'http://aug2020.archive.ensembl.org/',
  '100': 'http://apr2020.archive.ensembl.org/',
  '99': 'http://jan2020.archive.ensembl.org/',
  '98': 'http://sep2019.archive.ensembl.org/',
  '97': 'http://jul2019.archive.ensembl.org/',
  '96': 'http://apr2019.archive.ensembl.org/',
  '95': 'http://jan2019.archive.ensembl.org/',
  '94': 'http://oct2018.archive.ensembl.org/',
  '93': 'http://jul2018.archive.ensembl.org/',
  '92': 'http://apr2018.archive.ensembl.org/',
  '91': 'http://dec2017.archive.ensembl.org/',
  '90': 'http://aug2017.archive.ensembl.org/',
  '89': 'http://may2017.archive.ensembl.org/',
  '88': 'http://mar2017.archive.ensembl.org/',
  '87': 'http://dec2016.archive.ensembl.org/',
  '86': 'http://oct2016.archive.ensembl.org/',
  '80': 'http://may2015.archive.ensembl.org/',
  '77': 'http://oct2014.archive.ensembl.org/',
  '75': 'http://feb2014.archive.ensembl.org/',
  '54': 'http://may2009.archive.ensembl.org/'}

def test_ensembl_host(scplus_obj, host, species):
  dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
  annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
  annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
  annot['Chromosome'] = annot['Chromosome'].astype('str')
  filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
  annot = annot[~filter]
  annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
  gene_names_release = set(annot['Gene'].tolist())
  ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
  print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
  return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
  print(f'host: {version}')
  try:
    n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
  except:
    print('Host not reachable')

v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

