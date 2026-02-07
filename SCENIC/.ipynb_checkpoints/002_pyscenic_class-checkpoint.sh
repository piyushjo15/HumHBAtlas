#!/bin/bash

## this script generated adjacency list from pyscenic run for each class using each split
DIRF=$HOME/SCENICfiles
cd $HOME/SCENICOut/Class/$SCE
module load Miniconda3
source activate pyscn
echo "running pyscenic for $SCE"
#adjacency matrix
pyscenic grn -m genie3 ${SCE}_A.loom $DIRF/allTFs_hg38.txt -o adj${SCE}_A.csv --num_workers 6 --seed 10
#add correlation
pyscenic add_cor adj${SCE}_A.csv ${SCE}_A.loom --output cor${SCE}_A.tsv --mask_dropouts

pyscenic grn -m genie3 ${SCE}_B.loom $DIRF/allTFs_hg38.txt -o adj${SCE}_B.csv --num_workers 6 --seed 10
pyscenic add_cor adj${SCE}_B.csv ${SCE}_B.loom --output cor${SCE}_B.tsv --mask_dropouts
