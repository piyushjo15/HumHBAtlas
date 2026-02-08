import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle
##using generated models to predict MP class of CREs identified in fetal brain atlas
##local location
dataDir='/Deep/Outs/'
os.chdir(dataDir)

##load bed files, using augmented peaks
bed_files='FetalBrain/markerpeaks/'
chromsizes='extrafiles/chromSizes_canchrom.genome'
adata_ext = crested.import_beds(
    beds_folder=bed_files, chromsizes_file=chromsizes
)  
adata_ext
df=adata_ext.var
print(df.head())
df.to_csv('FetalBrain/MarkerCREsbed.txt',sep="\t")

import keras
genome_fasta=dataDir+'extrafiles/gencdH38p13r37.fa'

##load model
model_path = dataDir+"extrafiles/deephb_lstmv2_57MP_73.keras"

model = keras.models.load_model(
    model_path, compile=False
)  #
prediction = crested.tl.predict(adata_ext, model,genome=genome_fasta)

 np.savetxt('FetalBrain/pred_FB_deephb_lstmv2_57MP_73.txt', prediction, delimiter='\t', fmt='%.3f')
with open(os.path.join(out_Dir,'FetalBrain/pred_FB_deephb_lstmv2_57MP_73.pkl'), 'wb') as f:
        pickle.dump(prediction, f)

exit()
