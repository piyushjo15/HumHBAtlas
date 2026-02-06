#!/bin/bash

#PBS -N dedupbams
#PBS -l nodes=1:ppn=6
#PBS -l walltime=72:10:0
#PBS -l mem=60GB

DIR=$HOME/ATACana/Outs/CRI/
cd $DIR/dedup


module load SAMtools

#samtools based
samtools view -@ 6 -f 1 -F 1284 -q 20 -o ${SCE}.sorted.dedup.bam $DIR/bams/${SCE}.sorted.bam 
samtools index -@ 6 ${SCE}.sorted.dedup.bam
samtools flagstat ${SCE}.sorted.dedup.bam > ${SCE}.dedup.stat
