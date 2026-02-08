import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle


## gpu location
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
) 

#subset annData to remove row from MP00, the random background signature
#print(adata.obs.index)
adata=adata[adata.obs.index != 'MP00_CREs_aug']

# Standard train/val/test split
crested.pp.train_val_test_split(
    adata, strategy="region",random_state=145,val_size=0.11,test_size=0.0
)
print(adata.var["split"].value_counts())
## saving peaks that were assigned to test cohort
df=adata.var
df_test = df[df['split'] == "test"]
with open('deephb_testsplit_cnn_57MP.pkl', 'wb') as f:
    pickle.dump(df_test, f)

    
print("Configuring the model...")
# Datamodule
genome_fasta=dataDir+'extrafiles/gencdH38p13r37.fa'
datamodule = crested.tl.data.AnnDataModule(
    adata,
    genome=genome_fasta,
    batch_size=128,  
    max_stochastic_shift=0,  # This augemntation is not required as I am already augmenting
    always_reverse_complement=True, 
)


model_architecture = crested.tl.zoo.deeptopic_cnn(seq_len=500, 
                                                   num_classes=57,
                                                   filters=1024, ## default 1024
                                                   first_kernel_size=30, ##changed to 30
                                                   dense_out=1024, ## default 1024
                                                   pre_dense_do=0.5, #default
                                                   dense_do=0.5 #default
                                                   )

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
    run_name="CNN_57MP",
    logger=None,  
    seed=145
)

# train the model
trainer.fit(epochs=100) ## 100 is fine
met_test=trainer.test(return_metrics=True)
print('saving test metric as pickle..')
with open('met_deephb_cnn_57MP.pkl', 'wb') as f:
    pickle.dump(met_test, f)

exit()
