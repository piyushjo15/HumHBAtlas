#!/bin/bash

#BSUB -J homemotifana
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -W 2:00
#BSUB -n 8
#BSUB -R "rusage[mem=20GB]"


##script for running homer motif enrichment
#cd $HOME/MBsnANA/HBana/ATACana/Outs/TopicMP/TopicMP_CREs_57MP

#OUT=$HOME/MBsnANA/HBana/ATACana/Outs/TopicMP/TopicMP_CREs_57MP/motifenr
TMP=$HOME/MBsnANA/TMP

## all Hindbrain peaks
BG=$HOME/MBsnANA/HBana/ATACana/Outs/PeakCom/AllPeaks_robust.fil.sorted.bed

#cd $HOME/MBsnANA/HBana/ATACana/Outs/markerCREs/Class
#OUT=$HOME/MBsnANA/HBana/ATACana/Outs/markerCREs/Class/motifenr

cd $HOME/MBsnANA/HBana/ATACana/Outs/markerCREs/Pro_cluster
OUT=$HOME/MBsnANA/HBana/ATACana/Outs/markerCREs/Pro_cluster/motifenr
## combined peaks from all the MP
#BG=$HOME/MBsnANA/HBana/ATACana/Outs/TopicMP/ALLCREs_57MP.sort.mer.bed

#JAS=$HOME/MBsnANA/HBana/ATACana/Outs/extrafiles/JASPAR2024_root.motif
JAS=$HOME/MBsnANA/HBana/ATACana/Outs/extrafiles/allpat_lstmv2_57MP_73_220725.motif 

module load homer/4.11
echo "Processing motif analysis for Topic: $SCE"
##default
#findMotifsGenome.pl ${SCE}_CREs.bed hg38r $OUT/${SCE} -p 8 -size 500 -preparsedDir $TMP -fdr
## differential enriched comapred to all peaks
#findMotifsGenome.pl ${SCE}_CREs.bed hg38r ${OUT}2/${SCE} -p 8 -bg $BG -size 500 -preparsedDir $TMP -fdr

#findMotifsGenome.pl ${SCE}_CREs.bed hg38r \
#${OUT}_DH2/${SCE} -mknown $JAS -p 8 -size 500 -preparsedDir $TMP -fdr -nomotif ## no de novo ##no BG in DH2

#findMotifsGenome.pl ${SCE}_CREs.bed hg38r ${OUT}_DH/${SCE} -bg $BG \
#-mknown $JAS -p 8 -size 500 -preparsedDir $TMP -fdr   

findMotifsGenome.pl ${SCE}_CREs.bed hg38r ${OUT}_DH/${SCE} -bg $BG \
-mknown $JAS -p 8 -size 250 -preparsedDir $TMP -fdr   
#findMotifsGenome.pl ${SCE}_CREs.ext.v2.bed hg38r $OUT/${SCE}_extv2 -p 8 -size 500 -preparsedDir $TMP -fdr
# ###JASPAR
# findMotifsGenome.pl ${SCE}_CREs.sort.bed hg38r $OUT/${SCE}_JASPAR -mknown $JAS -bg $BG -p 8 -size 500 -preparsedDir $TMP -fdr 
# 
# #
# ###using annotatePeaks.pl to find which peaks have which mtoifs
# #annotatePeaks.pl beds/${SCE}.sorted.bed hg38r -m JASPAR2024.motif -nmotifs > $OUT/${SCE}_JASPAR2024motifenr.txt
# 
## -nmotifs = number of motif peaks
## -matrix <prefix> outputs a matrix couccurence 

# 
# ##script for running homer motif enrichment
# cd $HOME/MBsnANA/HBana/ATACana/Outs/markerCREs/
# 
# OUT=$HOME/MBsnANA/HBana/ATACana/Outs/markerCREs/motifenr
# TMP=$HOME/MBsnANA/TMP
# BG=$HOME/MBsnANA/HBana/ATACana/Outs/PeakCom/AllPeaks_robust.fil.sorted.bed
# JAS=$HOME/MBsnANA/HBana/ATACana/Outs/extrafiles/JASPAR2024.motif
# 
# module load homer/4.11
# echo "Processing motif analysis for Topic: $SCE"
# ###JASPAR
# findMotifsGenome.pl ${SCE}_CREs.bed hg38r $OUT/${SCE}_JASPAR -mknown $JAS -bg $BG -p 8 -size 500 -preparsedDir $TMP -fdr 
