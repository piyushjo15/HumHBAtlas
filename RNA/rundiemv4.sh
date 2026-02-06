#!/bin/bash

#PBS -N preprocessSNlibs
#PBS -l nodes=1:ppn=4
#PBS -l walltime=40:10:0
#PBS -l mem=50GB


cd $HOME/DATA/Diemana

module load R/4.2.2-foss-2022a
Rscript rundiemv4.R $DIEM ##HBlibsv1.txt
#Rscript rundiemv4b.R $DIEM ##for HBlibsv2.txt
Rscript runSoupx.R $DIEM

#Rscript runscrublet.R $DIEM ##run on command line
#Rscript runsceprov3.R $DIEM
#Rscript runQCplots.R $DIEM

