#!/bin/bash


##script for running homer motif enrichment
cd $HOME/ATACana/Outs/TopicMP/TopicMP_CREs_57MP

OUT=$HOME/ATACana/Outs/TopicMP/TopicMP_CREs_57MP/motifenr
TMP=$HOME/TMP
### This motif file obtained from syntaxes obtained from TFMoDisco analys
JAS=$HOME/ATACana/Outs/extrafiles/allpat_lstmv2_57MP_73_220725.motif 

## all Hindbrain peaks
BG=$HOME/MBsnANA/HBana/ATACana/Outs/PeakCom/AllPeaks_robust.fil.sorted.bed
module load homer/4.11
## FInding enrichment of identified motifs in each of the meta-regulatory programs
echo "Processing motif analysis for Topic: $SCE"
findMotifsGenome.pl ${SCE}_CREs.bed hg38r ${OUT}/${SCE} -bg $BG \
-mknown $JAS -p 8 -size 250 -preparsedDir $TMP -fdr 


## Ths analysis is for identifying motif enrichment in marker peaks
## associated with Progenitor clusters
cd $HOME/ATACana/Outs/markerCREs/Pro_cluster
OUT=$HOME/ATACana/Outs/markerCREs/Pro_cluster/motifenr
##marker peaks were identified per clusters within a class
echo "Processing motif analysis for cluster: $SCE"
findMotifsGenome.pl ${SCE}_CREs.bed hg38r ${OUT}/${SCE} -bg $BG \
-mknown $JAS -p 8 -size 250 -preparsedDir $TMP -fdr 