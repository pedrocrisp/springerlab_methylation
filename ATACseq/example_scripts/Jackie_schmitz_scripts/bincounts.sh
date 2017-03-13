#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/wenli_data/e_o/bincounts_o
#PBS -e /scratch.global/nosha003/wenli_data/e_o/bincounts_e
#PBS -N bincounts
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

cd /scratch.global/nosha003/schmitz/atacseq

## all bins
perl /home/springer/nosha003/wenli_data/scripts/bincounts.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}.100bpcounts.bed -o counts/${ID}.counts.txt

## macs overlapping bins --> only want bins in the intersect file (exclude zero counts)
perl /home/springer/nosha003/wenli_data/scripts/bincounts.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}_macs_bins_intersect.bed -o counts/${ID}_macs_bins_counts.txt

awk '$2 != 0' counts/${ID}_macs_bins_counts.txt > counts/${ID}_macs_bins_counts2.txt

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/bincounts.sh