#!/bin/bash

#PBS -N STARsolo
#PBS -l nodes=1:ppn=8
#PBS -l walltime=30:10:0
#PBS -l mem=60GB


## This script aligns dempltiplexed snRNA-seq reads generated 
## from for 10x single-cell RNAsequecncing. Reads were demultiplexed
## by the core facility using bcl2fastq

module load SAMtools/1.9-foss-2017a
module load STAR/2.7.9a

IDX=$HOME/index/starsologencdH38_100
OUT=$HOME/STARsolo
WL=$HOME/ann/3M-february-2018.txt
FQ=$HOME/fastqs/$SCE
cd $FQ
C1=`ls -m *${SCE}*_R1_001.fastq.gz | tr -d '\n' | tr -d ' '`
C2=`ls -m *${SCE}*_R2_001.fastq.gz | tr -d '\n' | tr -d ' '`

STAR --runMode alignReads --soloType CB_UMI_Simple \
--quantMode GeneCounts --runThreadN 8 \
--soloUMIlen 12 --soloCBwhitelist $WL \
--readFilesIn $C2 $C1 \
--genomeDir $IDX \
--twopassMode Basic \
--outFileNamePrefix $OUT/$SCE/$SCE \
--readFilesCommand zcat \
--soloFeatures Gene GeneFull \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--outSAMtype BAM Unsorted \
--soloCellFilter None \
--outSAMunmapped Within \
--outSAMmultNmax 1 \
--limitSjdbInsertNsj 1500000
