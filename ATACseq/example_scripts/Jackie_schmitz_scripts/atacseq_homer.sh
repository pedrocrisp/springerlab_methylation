#24Jan2017
#atacseq_homer.sh
#!/bin/bash

#HOMER peak calling
#PBS -l walltime=08:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_homer_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_homer_e
#PBS -N atacseq_homer
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

# install homer (http://homer.salk.edu/homer/introduction/install.html)
module load bedtools

#/home/springer/nosha003/software/homer/bin/makeTagDirectory /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_tag /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.bed -format bed
 
#/home/springer/nosha003/software/homer/bin/findPeaks /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_tag -region -size 200 -minDist 50 -o auto -tbp 0
                 
#/home/springer/nosha003/software/homer/bin/pos2bed.pl /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_tag/peaks.txt > /scratch.global/nosha003/schmitz/atacseq/homer/${ID}.tmp.bed

#bedtools sort -i /scratch.global/nosha003/schmitz/atacseq/homer/${ID}.tmp.bed | bedtools merge -i - > /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}.pks.bed


# uniq
/home/springer/nosha003/software/homer/bin/makeTagDirectory /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.bed -format bed
 
#/home/springer/nosha003/software/homer/bin/findPeaks /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag -region -size 200 -minDist 50 -o auto -tbp 0
/home/springer/nosha003/software/homer/bin/findPeaks /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag -region -size 200 -minDist 200 -o auto -tbp 0
                 
/home/springer/nosha003/software/homer/bin/pos2bed.pl /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq_tag/peaks.txt > /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq.tmp.bed

bedtools sort -i /scratch.global/nosha003/schmitz/atacseq/homer/${ID}_uniq.tmp.bed | bedtools merge -i - > /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_uniq.pks.bed

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_homer.sh