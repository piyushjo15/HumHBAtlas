#!/bin/bash

#PBS -N mergeBAM
#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:10:0
#PBS -l mem=80GB


module load SAMtools
module load BEDTools

DIR=$HOME/MBsnANA/HBana/ATACana/Outs/Superenhancer/

## loop version
##
# while read SCE
# do
#   echo "procesising Class $SCE.."
#   cd $DIR/$SCE
#   C1=`ls *_sp.bam`
#   samtools merge -@ 8 ${SCE}.merge.bam $C1
#   samtools sort -@ 8 ${SCE}.merge.bam -o ${SCE}.merge.sort.bam
#   samtools index -@ 8 ${SCE}.merge.sort.bam
#   cd $DIR
# done < $HOME/MBsnANA/HBana/ATACana/HBATAC_Scripts/$File

echo "procesising Class $SCE.."
cd $DIR/$SCE
C1=`ls *_sp.bam`
samtools merge -@ 8 ${SCE}.merge.bam $C1
samtools sort -@ 8 ${SCE}.merge.bam -o ${SCE}.merge.sort.bam
samtools index -@ 8 ${SCE}.merge.sort.bam
