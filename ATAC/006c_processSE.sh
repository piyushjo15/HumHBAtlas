#!/bin/bash
## procesed output from SE
cp ${line}_A_SuperStitched_REGION_TO_GENE.txt ${line}_A_super.bed
cp ${line}_B_SuperStitched_REGION_TO_GENE.txt ${line}_B_super.bed
sed -i '1d' ${line}_A_super.bed
sed -i '1d' ${line}_B_super.bed
## First obtain a bed file
awk -v OFS="\t" '{print $2, $3, $4, $7-$8}' ${line}_A_super.bed > ${line}_A_super.v2.bed
awk -v OFS="\t" '{print $2, $3, $4, $7-$8}' ${line}_B_super.bed > ${line}_B_super.v2.bed

# sort
while read line; do sort -k1,1 -k2,2n ${line}_A_super.v2.bed >${line}_A_super.v2.sort.bed; done < RNA_Class.txt
while read line; do sort -k1,1 -k2,2n ${line}_B_super.v2.bed >${line}_B_super.v2.sort.bed; done < RNA_Class.txt

# intsect
while read line; do bedtools intersect -a ${line}_A_super.v2.sort.bed -b ${line}_B_super.v2.sort.bed > ${line}_SE.bed; done < RNA_Class.txt