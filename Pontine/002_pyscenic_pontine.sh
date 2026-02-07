#!/bin/bash


## this script generated adjacency list from pyscenic run for each class using each split
DIRF=$HOME/SCENICfiles
cd $HOME//SCENICOut/PN
module load Miniconda3
source activate pyscn
echo "running pyscenic for pontine nuclei"
# ##using loom file generated from prepforSCENIC to get adjcancey matrcex
pyscenic grn -m genie3 PN.loom $DIRF/allTFs_hg38.txt -o adjPN.csv --num_workers 6 --seed 10

echo "generating correlation matrix for pontine nucle"

#add correlation
pyscenic add_cor adjPN.csv PN.loom --output corPN.tsv --mask_dropouts
