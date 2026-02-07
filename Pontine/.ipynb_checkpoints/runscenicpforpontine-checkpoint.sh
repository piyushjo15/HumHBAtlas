#!/bin/bash

#BSUB -J runscenicpontine
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 52:00
#BSUB -n 8
#BSUB -R "rusage[mem=130GB]"

cd ~/MBsnANA/HBana/SCENICana/SCENIC_Scripts
module load miniconda/4.9.2
source activate scenicplus

#python 002b_runTopicModeling_pontine.py
#python 003_runPycisTarget_pontine.py
#python 005a_scplusobj_pontine.py
#python 006a_runSCENIC_pontine.py
python 006a_runSCENIC_pontinev2.py