#!/bin/bash

#MSUB -A b1042
#MSUB -q genomicslong
#MSUB -l walltime=72:00:00
#MSUB -M behram.radmanesh@northwestern.edu
#MSUB -l nodes=1:ppn=20
#MSUB -e errlog

cd /projects/p30459/guidescan/gencodev22_run/

source activate guidescan

guidescan_guidequery -b ../cas9_hg38_all_guides.bam --batch gencode_v22_lncRNA_transcripts.gtf --target flanking --flankdistance 1000 --sort offtargets -n 7 --annot gencodev22_all_annotations.bed --blat ~/miniconda3/envs/guidescan/bin/blat -o guidescan_out
