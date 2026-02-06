#!/bin/bash

#PBS -N STARsoloMSN
#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:10:0
#PBS -l mem=45GB
module load SAMtools/1.9-foss-2017a
export PATH=$PATH:$HOME/STAR/source

IDX=$HOME/HB/index/starsologencdH38_100
OUT=$HOME/DATA/STARsolo
WL=$HOME/DATA/CR_ARC/lib/python/cellranger/barcodes/737K-arc-v1.txt
FQ=$HOME/HB/HMoPo/fastqs/${SCE}
cd $FQ
C1=`ls -m *${SCE}*_R1_001.fastq.gz | tr -d '\n' | tr -d ' '`
C2=`ls -m *${SCE}*_R2_001.fastq.gz | tr -d '\n' | tr -d ' '`

STAR --runMode alignReads --soloType CB_UMI_Simple \
--quantMode GeneCounts --runThreadN 8 \
--soloUMIlen 12 --soloCBwhitelist $WL \
--readFilesIn $C2 $C1 \
--genomeDir $IDX \
--twopassMode Basic \
--outFileNamePrefix $OUT/${SCE}/${SCE} \
--readFilesCommand zcat \
--soloFeatures Gene GeneFull \
--soloUMIfiltering MultiGeneUMI \
--soloCBmatchWLtype 1MM_multi_pseudocounts \
--outSAMtype BAM Unsorted \
--soloCellFilter None \
--outSAMunmapped Within \
--outSAMmultNmax 1 \
--limitSjdbInsertNsj 1500000
