#!/bin/bash

#PBS -N macs3run
#PBS -l nodes=1:ppn=7
#PBS -l walltime=1:10:0
#PBS -l mem=10GB

module load Bowtie2
module load SAMtools
module load Miniconda3


DIR=$HOME/ATACana/Outs/CRI/
TMP=$HOME/TMP

source ~/conda_env/packages/macs3/bin/activate
cd $IN
echo "calling narrow peaks for $SCE"


##extract over lapping CREs from macs3 called peaks
echo "post processing after MACS3 peak call...$SCE."
AllPeak=$HOME/ATACana/Outs/PeakCom/AllPeaks_robust.fil.sorted.bed
module load BEDTools
cd ${SCE}_A
#sort
sort -k1,1 -k2,2n ${SCE}_A_peaks.narrowPeak > ${SCE}_A.sort.bed
##intersect
bedtools intersect -wa -a $AllPeak -b ${SCE}_A.sort.bed > ${SCE}_A_peaks4SE.bed

cd ..
cd ${SCE}_B
#sort
sort -k1,1 -k2,2n ${SCE}_B_peaks.narrowPeak > ${SCE}_B.sort.bed
##intersect
bedtools intersect -wa -a $AllPeak -b ${SCE}_B.sort.bed > ${SCE}_B_peaks4SE.bed



