import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle

##read id (e.g. MP10)
MP=sys.argv[1]
## gpu location
dataDir='/Deep/'
dataDir='/Deep/Outs/'

os.chdir(dataDir)


##load bed files, using non-augmentnted peaks
## running test for one
bed_files=dataDir+'TopicMP_CREs/'
chromsizes=dataDir+'extrafiles/chromSizes_canchrom.genome'
fn=MP+'_CREs'
adata_ext = crested.import_beds(
    beds_folder=bed_files, chromsizes_file=chromsizes,
    classes_subset=[fn]
)  
adata_ext

##load model
model_path = dataDir+'Outs/deephb_lstmv2_57MP_73.keras"

import keras

model = keras.models.load_model(
    model_path, compile=True
)  # 
##find index of metaprogram
with open(os.path.join(dataDir,'extrafiles/TopicMPs.txt'), "r") as file:
    MPs = file.read().splitlines()
print(MPs)
idx=MPs.index(MP)

##calculate contribution scores
genome_fasta=dataDir+'extrafiles/gencdH38p13r37.fa'

scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
    adata_ext,
    target_idx=idx,  # None (=all classes), list of target indices, or empty list (='combined' class)
    model=model,
    method="expected_integrated_grad",  # default. Other options: "integrated_grad", "mutagenesis"
    genome=genome_fasta,
    all_class_names =MPs,
    transpose=True,
    output_dir=os.path.join(out_Dir,'cs_lstmv2_57MP_73')
)

##run tfmodisco

crested.tl.modisco.tfmodisco(
    contrib_dir =os.path.join(out_Dir,'cs_lstmv2_57MP_73'),
    class_names=[MP],
    output_dir=os.path.join(out_Dir,'mc_lstmv2_57MP_73'),
    sliding_window_size =15,
    min_metacluster_size =100)

print("Finished running TF-modisco for metaprogram "+MP)