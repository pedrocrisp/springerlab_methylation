#13Feb17
#atacseq_filter_fpkm.sh
#!/bin/bash

#filter peaks
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=20gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_filter_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_filter_e
#PBS -N atacseq_filter
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"
READS="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 5)"

module load bedtools

bedtools intersect -a /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_uniq.pks.bed -b /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.bed -wo -bed > /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_reads.pks.bed

perl /home/springer/nosha003/schmitz/scripts/homer_peaks_counts.pl -i /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_reads.pks.bed -I /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_reads.pks.bed -o /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_counts.pks.txt 

perl /home/springer/nosha003/schmitz/scripts/atac_homer_fpkm_filter.pl -n ${READS} -i /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_counts.pks.txt -o /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_peaks_fpkm.bed

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_filter2.sh
# qsub -t 1-4 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_combinereps_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_filter2.sh