#!/bin/bash

#PBS -N run_rose
#PBS -l nodes=1:ppn=4
#PBS -l walltime=48:10:0
#PBS -l mem=130GB

export PATH=$PATH:$HOME/ATACana/ROSE/bin
cd $HOME/ATACana/ROSE/bin

module load BEDTools
module load SAMtools
module load R


IN=$HOME/ATACana/Outs/CRI/
OUT=$HOME/ATACana/Outs/CRI/rose
BGD=$HOME/ATACana/Outs/CRI/BGD/BGD.merge.sort.bam

BAM=$IN/${SCE}/${SCE}.merge.sort.bam

PEAK=$IN/${SCE}/${SCE}_peaks4SE.bed
GFF=$IN/${SCE}/${SCE}.gff

echo "generating GFF from bed"
awk '{OFS="\t"; print $1, $4, ".", $2, $3, ".",".",".", $4}' ${PEAK} > ${GFF}

echo "Running ROSE for ..."
echo $SCE
./ROSE_main.py -g HG38 -i $GFF -r $BAM -o ${OUT} -t 2000 -c $BGD

