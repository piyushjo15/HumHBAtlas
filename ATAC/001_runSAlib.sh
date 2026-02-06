#!/bin/bash

#PBS -N runSAalign
#PBS -l nodes=1:ppn=8
#PBS -l walltime=40:10:0
#PBS -l mem=200GB

## aligning ATAC reads using cellranger-atac
export PATH=$PATH:$HOME/CRatac
cd $HOME/snATACana/

cellranger-atac count --id=${SCE} \
		--reference=$HOME/index/CR_ARC_HG38 \
		--fastqs=$HOME/fastqs/${SCE} \
		--sample=${SCE} 
