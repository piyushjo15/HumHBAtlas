#!/bin/bash

#PBS -N splitbam
#PBS -l nodes=1:ppn=6
#PBS -l walltime=20:10:0
#PBS -l mem=2GB


##using new function to subset bams using 10x provided app
#https://github.com/10XGenomics/subset-bam?tab=readme-ov-file
module load SAMtools
module load BEDTools

echo "processing Class n split $SCE"
subset=$HOME/conda_env/packages/subset-bam_linux
DIR=$HOME/ATACana/Outs/CRI/

cd $DIR/$SCE

while read line
do
  echo "processing for Sample $line"
  #cd $DIR
  mkdir -p $HOME/TMP/${SCE}_${line}
  export TMPDIR=$HOME/TMP/${SCE}_${line} #define tmp dir
  $subset --bam $DIR/dedup/${line}.sorted.dedup.bam --cell-barcodes ${line}_Index.csv --cores 6 --out-bam ${line}_sp.bam
done < Samples.txt

# ## this for background
# echo "procesising sample $SCE"
# DIR=$HOME/ATACana/Outs/CRI/
# DIRB=$HOME/ATACana/Outs/CRI/BGD/ # for background
# cd $DIRB
# #make OUT dirs
# mkdir -p $HOME/TMP/$SCE
# export TMPDIR=$HOME/TMP/$SCE #define tmp dir
# 
# $subset --bam $DIR/dedup/${SCE}.sorted.dedup.bam \
# --cell-barcodes ${SCE}_Index.csv --cores 6 --out-bam ${SCE}_bgd.bam
# ## i had to add "-1" to each of the barcodes using while loop
# #sed 's/$/-1/' ${SCE}_Index.csv > $outputname