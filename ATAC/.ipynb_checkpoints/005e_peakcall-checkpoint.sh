#!/bin/bash

#PBS -N peakcall
#PBS -l nodes=1:ppn=6
#PBS -l walltime=4:10:0
#PBS -l mem=60GB

##this script first utilzies deeptools to shift ATAc reads and then MACS2 to call peaks

DIR=$HOME/ATACana/Outs/CRI/
BL=$HOME/ann/hg38-blacklist.v2_noncanchr.bed
PD=$HOME/ATACana/Outs/CRI/broadpeaks

cd $DIR


module load Miniconda3
module load SAMtools
module load BEDTools
module load MACS2


## parallel version
cd $SCE

#echo "removing unwanted bam file  .."
rm ${SCE}_merge.bam

##running this without peak call
echo "converting bam file to BED .."

macs2 randsample -i ${SCE}_merge.sort.bam -f BAMPE -p 100 --seed 1 -o ${SCE}_merge.bed

echo "calling peaks .."

###peak call

macs2 callpeak -f BEDPE --nomodel --nolambda --shift -75 --extsize 150 --seed 1 \
-g hs --broad --keep-dup all \
-q 0.001 -n ${SCE} -t ${SCE}_merge.bed --outdir $PD 2> $PD/${SCE}_macs2.log

cd $PD
echo "sort peak file and removing blacklisted regions .."
sort -k1,1 -k2,2n ${SCE}_peaks.broadPeak > ${SCE}_bp.sort.bed
##remove blacklisted region and unwanted chromosomes
bedtools subtract -a ${SCE}_bp.sort.bed -b $BL > ${SCE}_bp.fil.sort.bed
