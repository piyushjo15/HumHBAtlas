## This script runs topic modeling on cisTopic object

import os
import pandas
import numpy
import warnings
import pycisTopic
import sys
import pickle
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *

#Define dir

projDir = '/home/SCENICOut/PN/Cistopic'
os.chdir(projDir)
tmpDir = '/home/TMP/'

cistopic_obj = pickle.load(open(os.path.join(projDir, 'PN_cistopic_obj.pkl'), 'rb'))

print("running topic modelling for 50 topics")
models=run_cgs_models(cistopic_obj,
                        n_topics=[50],
                        n_cpu=4,
                        n_iter=500,
                        random_state=555,
                        alpha=50,
                        alpha_by_topic=True,
                        eta=0.1,
                        eta_by_topic=False,
                        save_path=None,
                        _temp_dir = tmpDir)

if not os.path.exists(os.path.join(projDir, 'models')):
  os.makedirs(os.path.join(projDir, 'models'))

pickle.dump(models, open(os.path.join(projDir, 'models/ATAC_models_500_iter_LDA.pkl'), 'wb'))
##loading models
#models = pickle.load(open(os.path.join(projDir, 'models/ATAC_models_500_iter_LDA.pkl'), 'rb'))

print("evluating and adding topic models")
model=evaluate_models(models,
                      select_model=50,
                      return_model=True,
                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                      plot_metrics=False)

# Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,open(os.path.join(projDir, 'PN_cistopic_obj.pkl'), 'wb'))

#First we will binarize the topics using the otsu method and by taking the top 3k regions per topic.
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

#Run without imputation
#Run PycisTopic without imputation
nonimputed_obj = CistopicImputedFeatures(imputed_acc=cistopic_obj.fragment_matrix,cell_names=cistopic_obj.cell_names,feature_names=cistopic_obj.region_names,project="cistopic_not_imputed")

# remove zero rows as done by the impute_accessibility function
zero_rows = (nonimputed_obj.mtx.sum(axis=1) == 0)
nonimputed_obj.mtx = nonimputed_obj.mtx[~np.array(zero_rows).flatten(), :]
nonimputed_obj.feature_names = list(np.array(nonimputed_obj.feature_names)[~np.array(zero_rows).flatten()])


#Save Results
if not os.path.exists(os.path.join(projDir, 'candidate_enhancers')):
  os.makedirs(os.path.join(projDir, 'candidate_enhancers'))

pickle.dump(region_bin_topics_otsu, open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))

#save cistopic object
pickle.dump(cistopic_obj, open(os.path.join(projDir, 'PN_cistopic_obj.pkl'), 'wb'))

print("Finished 002_selectModel for pontine")



exit()


