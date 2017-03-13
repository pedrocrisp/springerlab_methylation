#Feb_7_2017
#atacseq_combinereps_homer
#!/bin/bash

#Bedgraph 
#PBS -l walltime=08:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/homer_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/homer_e
#PBS -N homer
#PBS -V
#PBS -r n

ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST})"

module load bedtools

cd /scratch.global/nosha003/schmitz/atacseq/align
/home/springer/nosha003/software/homer/bin/makeTagDirectory /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.bed -format bed
 
/home/springer/nosha003/software/homer/bin/findPeaks /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag -region -size 200 -minDist 50 -o auto -tbp 0
                 
/home/springer/nosha003/software/homer/bin/pos2bed.pl /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag/peaks.txt > /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq.tmp.bed

bedtools sort -i /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq.tmp.bed | bedtools merge -i - > /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_uniq.pks.bed

# qsub -t 1-4 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_combinereps_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_combinereps_homer.sh