import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle

## gpu location
##load bed files
dataDir='/dkfz/cluster/gpu/data/OE0290/p541i/'
out_Dir='/dkfz/cluster/gpu/checkpoints/OE0290/p541i/'
os.chdir(out_Dir)
## this script predicts the class os unaugmented Topic MP CREs
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

##load model
#model_path = dataDir+"extrafiles/deephb_cnn_99.keras"
#model_path = dataDir+"extrafiles/deephb_lstmv3_57MP_53.keras"
model_path = dataDir+"extrafiles/deephb_lstmv2_57MP_73.keras"


model = keras.models.load_model(
    model_path, compile=False
)  #

prediction = crested.tl.predict(adata_ext, model,genome=genome_fasta)

# ##multiple models
# model_path1 = dataDir+"extrafiles/deephb_lstmv2_05.keras"
# model_path2 = dataDir+"extrafiles/deephb_lstmv2_51.keras"
# model_path3 = dataDir+"extrafiles/deephb_lstmv2_95.keras"
# 
# model1 = keras.models.load_model(
#     model_path1, compile=False
# )  # 
# model2 = keras.models.load_model(
#     model_path2, compile=False
# )  #
# model3 = keras.models.load_model(
#     model_path3, compile=False
# )  #
# model_list=[model1,model2,model3]

#prediction = crested.tl.predict(adata_ext, model_list,genome=genome_fasta)


## save prediction

print("saving predictions..")
#with open(os.path.join(out_Dir,'pred_HBCREs_deephb_cnn_57MP_99.pkl'), 'wb') as f:
#      pickle.dump(prediction, f)

with open(os.path.join(out_Dir,'pred_HBCREs_deephb_lstmv2_57MP_73.pkl'), 'wb') as f:
      pickle.dump(prediction, f)

#np.savetxt('pred_HBCREs_deephb_lstmv2_95.txt', prediction, delimiter='\t', fmt='%.3f')
# with open(os.path.join(out_Dir,'pred_HBCREs_deephb_lstmv2_05.pkl'), 'wb') as f:
#        pickle.dump(prediction, f)
# with open(os.path.join(out_Dir,'pred_HBCREs_deephb_lstmv2_055195mix.pkl'), 'wb') as f:
#        pickle.dump(prediction, f)

#np.savetxt('pred_HBCREsext_deephb_lstmv2_95.txt', prediction, delimiter='\t', fmt='%.3f')
#with open(os.path.join(out_Dir,'pred_HBCREsext_deephb_lstmv2_95.pkl'), 'wb') as f:
#       pickle.dump(prediction, f)


exit()
