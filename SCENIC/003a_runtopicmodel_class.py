
import os
import pandas as pd
import numpy as np
import warnings
import pycisTopic
import pickle
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *

## creating cisTopic object for Pyscistarget modeling
#Define dir
Cls = sys.argv[1]
Splt=sys.argv[2] ## 'A' or 'B'
print("Starting select model for Class: " + Cls+" and split: "+Splt)

DIROUT='/home/SCENICOut/Class/'
projDir = os.path.join(DIROUT,Cls,'SCENIC',Cls+'_'+Splt)
os.chdir(projDir)
tmpDir='/home/TMP/'

####### MODEL SELECTION ##########
cistopic_obj = pickle.load(open(os.path.join(dataDir, Cls+'_'+Splt+'_cistopic_obj.pkl'), 'rb'))

projDir2 = projDir+"/models"
if not os.path.exists(projDir2):
  os.makedirs(projDir2)
print("Running models:")
## creat topic models
models=run_cgs_models(cistopic_obj,
                      n_topics=[50],
                      n_cpu=7,
                      n_iter=500,
                      random_state=555,
                      alpha=50,
                      alpha_by_topic=True,
                      eta=0.1,
                      eta_by_topic=False,_temp_dir = tmpDir,
                      save_path=None)
pickle.dump(models,open(os.path.join(projDir2, 'ATAC_models_500_iter_LDA.pkl'), 'wb'))

model=evaluate_models(models,
                      select_model=50,
                      return_model=True,
                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                      plot_metrics=False)

# Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)

# Save
pickle.dump(cistopic_obj,open(os.path.join(projDir, Cls+'_'+Splt+'_cistopic_obj.pkl'), 'wb'))


#First we will binarize the topics using the otsu method and by taking the top 2k regions per topic.
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top2k = binarize_topics(cistopic_obj, method='ntop', ntop = 2000)

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
pickle.dump(region_bin_topics_top2k, open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_top2k.pkl'), 'wb'))

#save cistopic object
pickle.dump(cistopic_obj,open(os.path.join(projDir, Cls+'_'+Splt+'_cistopic_obj.pkl'), 'wb'))

print("Finished 002_selectModel for:" + Cls+'_'+Splt)
exit()
