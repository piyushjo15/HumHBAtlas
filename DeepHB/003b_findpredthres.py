import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle

## gpu location
##load bed files
dataDir='/home/DeepHB/Outs'
os.chdir(dataDir)

## this script predicts the class of unaugmented Topic MP CREs
## I also use the output for auROC/PR values
## I also use extended sets of MP CREs for ROC/PR values
##load bed files, using augmented peaks
bed_files=dataDir+'TopicMP_CREs/'

chromsizes=dataDir+'extrafiles/chromSizes_canchrom.genome'
adata_ext = crested.import_beds(
    beds_folder=bed_files, chromsizes_file=chromsizes
)  # the regions file is optional for import_beds
adata_ext

# gt=adata_ext.X.T #ground turth
# print(gt.shape)

#with open(os.path.join(dataDir,'extrafiles/TopicMPs.txt'), "r") as file:
#    labels = file.read().splitlines()
#print(labels)


import keras
genome_fasta=dataDir+'extrafiles/gencdH38p13r37.fa'

##load model, the last best epoc was selected
#model_path = dataDir+"Outs/deephb_cnn_57MP_99.keras"
#model_path = dataDir+"Outs/deephb_lstmv1_57MP_53.keras"
model_path = dataDir+"Outs/deephb_lstmv2_57MP_73.keras"

model = keras.models.load_model(
    model_path, compile=False
)  #

prediction = crested.tl.predict(adata_ext, model,genome=genome_fasta)

## save prediction

print("saving predictions..")
#with open(os.path.join(out_Dir,'pred_HBCREs_deephb_cnn_57MP_99.pkl'), 'wb') as f:
#      pickle.dump(prediction, f)
with open(os.path.join(out_Dir,'pred_HBCREs_deephb_lstmv1_57MP_53.pkl'), 'wb') as f:
      pickle.dump(prediction, f)
with open(os.path.join(out_Dir,'pred_HBCREs_deephb_lstmv2_57MP_73.pkl'), 'wb') as f:
      pickle.dump(prediction, f)


exit()
