import os
import sys
import random
import numpy as np
import pandas as pd
import crested
import pickle

dataDir='/Deep/'
os.chdir(dataDir)

with open(os.path.join(dataDir,'TopicMPs.txt'), "r") as file:
    labels = file.read().splitlines()
print(labels)

modDir=dataDir+"mc_lstmv2_57MP_73/"

# First we obtain the resulting modisco files per class
matched_files = crested.tl.modisco.match_h5_files_to_classes(
    contribution_dir=modDir, classes=labels
)

# Then we cluster matching patterns, and define a pattern matrix [#classes, #patterns] describing their importance
all_patterns = crested.tl.modisco.process_patterns(
    matched_files,
    sim_threshold=4.25,  # The similarity threshold used for matching patterns. We take the -log10(pval), pval obtained through TOMTOM matching from tangermeme
    trim_ic_threshold=0.05,  # Information content (IC) threshold on which to trim patterns
    discard_ic_threshold=0.2,  # IC threshold used for discarding single instance patterns
    verbose=True,  # Useful for doing sanity checks on matching patterns
)
pattern_matrix = crested.tl.modisco.create_pattern_matrix(
    classes=labels,
    all_patterns=all_patterns,
    normalize=False,
    pattern_parameter="seqlet_count_log",
)
pattern_matrix.shape

print('saving results..')
print('saving results..')
with open('allpat_lstmv2_57MP_73.pkl', 'wb') as f:
    pickle.dump(all_patterns, f)


with open('pm_lstmv2_57MP_73.pkl', 'wb') as f:
    pickle.dump(pattern_matrix, f)
    
exit()