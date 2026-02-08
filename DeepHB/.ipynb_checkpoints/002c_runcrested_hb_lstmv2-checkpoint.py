import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle

## LSTM v2 is CNN+RNN hybrid model with JASPAR motif initialization
##load bed files
dataDir='/home/DeepHB/'
out_Dir='/home/DeepHB/Outs'
os.chdir(out_Dir)


##load bed files, using augmented peaks
## The peaks are augments by 50 bp windown upto 100 bp
bed_files=dataDir+'TopicMP_CREs_aug/'
chromsizes=dataDir+'extrafiles/chromSizes_canchrom.genome'
adata = crested.import_beds(
    beds_folder=bed_files, chromsizes_file=chromsizes
)  # the regions file is optional for import_beds
#adata

#subset annData to remove row from MP00, the random background signature
#print(adata.obs.index)
adata=adata[adata.obs.index != 'MP00_CREs_aug']

# Standard train/val/test split
crested.pp.train_val_test_split(
    adata, strategy="region",random_state=145
)
print(adata.var["split"].value_counts())

print("Configuring the model...")
# Datamodule
genome_fasta=dataDir+'extrafiles/gencdH38p13r37.fa'
datamodule = crested.tl.data.AnnDataModule(
    adata,
    genome=genome_fasta,
    batch_size=128,  
    max_stochastic_shift=0,  
    always_reverse_complement=True,
)


motif_dict=dataDir+"extrafiles/JASPAR2024_cluster_motifs_dict.pkl"## path to JASPAR 2024 PWMs

model_architecture = crested.tl.zoo.deeptopic_lstm(seq_len=500, 
                                                   num_classes=57,
                                                   filters=1024,
                                                   first_kernel_size=30, ## size of biggest motif is 29
                                                   max_pool_size=16,
                                                   max_pool_stride=8, 
                                                   dense_out=1024,
                                                   lstm_out=512, 
                                                   lstm_do=0.3, 
                                                   pre_dense_do=0.3, 
                                                   dense_do=0.4, 
                                                   motifs_path=motif_dict)

# Config: we will use the default topic classification config 
# (binary cross entropy loss and AUC/ROC metrics)
# This is appropriate
config = crested.tl.default_configs("topic_classification")
print(config)


print("Running training model...")

# setup the trainer
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=config,
    project_name="DeepHB", 
    run_name="LSTMv2_57MP", ## with motif intialization
    logger=None,  
    seed=145
)

# train the model
trainer.fit(epochs=100) ### only 100 epochs are fine

met_test=trainer.test(return_metrics=True)
print('saving test metric as pickle..')
with open('met_deephb_lstmv2_57MP.pkl', 'wb') as f:
    pickle.dump(met_test, f)

exit()
