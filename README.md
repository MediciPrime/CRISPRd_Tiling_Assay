## Analysis Overview
The point of this analysis is to ask GuideScan to find pgRNAs for all lncRNA
transcripts for Gencode v22 (GRCh38). This was done because the CRISPRa paper
used this version of GENCDOE and Jennie wants to utilize their activation list.
It also seems that the activation assay was performed on lncRNAs transcripts
because the activation approach has transcript names. Thus we will only use
transcript and disregard the genes that show up on the GTF file.

### Pipeline
1. Take Gencode v22 lncRNAs and extract transcript information.
   + `grep -P 'transcript\t' gencode_v22_lncRNA.gtf > gencode_v22_lncRNA_transcript.gtf`
2. GuideScan also needs a 'bed' annotation to prevent guide creation
   + Obtain Gencodev22 Primary Annotations and remove all mentions of lncRNAs
     + `sed -i.bak '/lncRNA/d' ./gencodev22_all_annotations.gtf`
     + `sed -i.bak '/lincRNA/d' ./gencodev22_all_annotations.gtf`
3. Convert the modified GTF file into a bed file for use with GuideScan
   + `awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencodev22_all_annotations.gtf | gtf2bed - > gencodev22_all_annotations.bed`
4. Modify the annotation bed file and extend UTR windows by 5kb to the left or right depending on strand orientation
   + The 'modify_anno' function is located in the 'pgRNA_creator.py' file, used to add 5kb windows
5. Run GuideScan using the following command:
   + `guidescan_guidequery -b ../cas9_hg38_all_guides.bam --batch gencode_v22_lncRNA_transcripts.gtf --target flanking --flankdistance 1000 --sort offtargets -n 7 --annot gencodev22_all_annotations.bed --blat ~/miniconda3/envs/guidescan/bin/blat -o guidescan_out`
6. Use the method 'pgRNA_clean' from 'pgRNA_creator.py' to create the pgRNA list ('pgRNA.txt') file.

### Initial Checks
There was a fear that a good portion of lncRNAs may intersect with other lncRNAs.
This fear was confirmed when the 'idOverlap' function found almost 31% of all
9,999 lncRNA targets were overlapping. The main problem with this would be the
possibility that the deletion of one lncRNA may lead to the deletion of a
neighboring lncRNA. However since we are also performing an activation CRISPR
assay, a check can be placed by using both the activation and deletion assays.
