#Jan_28_2017
#bincounts

#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/bincounts_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/bincounts_e
#PBS -N bincounts
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

cd /scratch.global/nosha003/schmitz/atacseq

## all bins
#perl /home/springer/nosha003/schmitz/scripts/bincounts_bed.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}.bowtie.100bpcounts.noheader.bed -o counts/${ID}.bowtie.counts.txt
#perl /home/springer/nosha003/schmitz/scripts/bincounts_bed.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I /home/springer/nosha003/schmitz/counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed -o counts/${ID}_uniq.bowtie.counts.txt

## macs overlapping bins --> only want bins in the intersect file (exclude zero counts)
#perl /home/springer/nosha003/schmitz/scripts/bincounts_bed.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}_bowtie_macs_bins.bed -o counts/${ID}_bowtie_macs_bins_counts.txt
#perl /home/springer/nosha003/schmitz/scripts/bincounts_bed.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}_uniq_bowtie_macs_bins.bed -o counts/${ID}_uniq_bowtie_macs_bins_counts.txt

## homer overlapping bins --> only want bins in the intersect file (exclude zero counts)
#perl /home/springer/nosha003/schmitz/scripts/bincounts_bed.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}_bowtie_homer_bins.bed -o counts/${ID}_bowtie_homer_bins_counts.txt
perl /home/springer/nosha003/schmitz/scripts/bincounts_bed_homer.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}_uniq_bowtie_homer_bins.bed -o counts/${ID}_uniq_bowtie_homer_bins_counts.txt



#awk '$2 != 0' counts/${ID}_bowtie_macs_bins_counts.txt > counts/${ID}_bowtie_macs_bins_counts2.txt
#awk '$2 != 0' counts/${ID}_bowtie_homer_bins_counts.txt > counts/${ID}_bowtie_homer_bins_counts2.txt
#awk '$2 != 0' counts/${ID}_uniq_bowtie_macs_bins_counts.txt > counts/${ID}_uniq_bowtie_macs_bins_counts2.txt
awk '$2 != 0' counts/${ID}_uniq_bowtie_homer_bins_counts.txt > counts/${ID}_uniq_bowtie_homer_bins_counts2.txt

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/bincounts_bowtie.sh