## This script runs extract ntop regions per topics for the array of ntopics I ran
## it is run in the conda environment scenicplus
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
DIR='/home/ATACana/Outs/TopicMP/'
Clss = sys.argv[1] #Class
Split="B" ## "A" or "B"
Clss = Clss+"_"+Split
projDir = os.path.join(DIR,Clss)
os.chdir(projDir)
tmpDir='/home/TMP/'

cistopic_obj = pickle.load(open(os.path.join(projDir, Clss+'_cistopic_obj.pkl'), 'rb'))

projDir2 = projDir+"/models"
if not os.path.exists(projDir2):
  os.makedirs(projDir2)
print("Running models:")
## creat topic models
models=run_cgs_models(cistopic_obj,
                      n_topics=[10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75],
                      n_cpu=7,
                      n_iter=500,
                      random_state=555,
                      alpha=50,
                      alpha_by_topic=True,
                      eta=0.1,
                      eta_by_topic=False,_temp_dir = tmpDir,
                      save_path=None)

pickle.dump(models,open(os.path.join(projDir2, 'ATAC_models_500_iter_LDA.pkl'), 'wb'))

my_list = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]

##load model and cistopic object
#models = pickle.load(open(os.path.join(projDir2, 'ATAC_models_500_iter_LDA.pkl'), 'rb'))

# Iterate through the list and print each element
if not os.path.exists(os.path.join(projDir, 'topics')):
  os.makedirs(os.path.join(projDir, 'topics'))

for item in my_list:
  print("Adding model for topic rank:")
  print(item)
  model=evaluate_models(models,select_model=item,return_model=True,metrics=['Minmo_2011'],plot_metrics=False)
  # Add model to cisTopicObject
  cistopic_obj.add_LDA_model(model)
  top2k = binarize_topics(cistopic_obj, method='ntop', ntop = 2000)
  #to save another way, just in case
  #pickle.dump(top2k, open(os.path.join(projDir, 'topics/topic',item,'_regions.pkl'), 'wb'))
  ## convert above list to dataframe and save
  top2k_df = pd.concat([pd.DataFrame({'DataFrame_Name':k, 'Row_Index':df.index, 'Value':df[k]}) for k, df in top2k.items()], ignore_index=True)
  top2k_df.to_csv(os.path.join(projDir, 'topics/topic'+str(item)+'_regions.txt'), sep='\t', index=False)
  


print("Finished extracting ntopics and regions for: " + Clss)


exit()


