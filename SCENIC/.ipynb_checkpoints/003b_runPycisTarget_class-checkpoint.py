
import os
import pandas
import numpy
import warnings
import pycistarget
import pickle

## Running Pycistarget on per class and split topics
#Define Project dir and temp dir
import sys

#Define dir
Cls = sys.argv[1]
Splt=sys.argv[2] ## 'A' or 'B'
print("Running PycisTarget for Class: " + Cls+" and split: "+Splt)

DIROUT='/home/SCENICOut/Class/'
projDir = os.path.join(DIROUT,Cls,'SCENIC',Cls+'_'+Splt)

tmpDir='/home/TMP/'


#load candidate enhancer regions identified in previous step, didn't use DARs
region_bin_topics_otsu = pickle.load(open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top2k = pickle.load(open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_top2k.pkl'), 'rb'))


#convert to pyranges objects
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for topic in region_bin_topics_top2k.keys():
    regions = region_bin_topics_top2k[topic].index[region_bin_topics_top2k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))


for key in region_sets.keys():
  print(f'{key}: {region_sets[key].keys()}')

#Define rankings, score and motif annotation database
db_fpath ='/home//SCENICfiles/'
motif_annot_fpath ='/home/SCENICfiles/'

rankings_db = os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')


from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
  region_sets = region_sets,
  species = 'homo_sapiens',
  save_path = projDir,
  ctx_db_path = rankings_db,
  dem_db_path = scores_db,
  path_to_motif_annotations = motif_annotation,
  run_without_promoters = True,
  n_cpu = 1,
  biomart_host = 'http://www.ensembl.org',
  _temp_dir = os.path.join(tmpDir,'ray_spill'),
  annotation_version = 'v10nr_clust',
)


print("Finished pycisTarget for for Class: " + Cls+" and split: "+Splt)

exit()

