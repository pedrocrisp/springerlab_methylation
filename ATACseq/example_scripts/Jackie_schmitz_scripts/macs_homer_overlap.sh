#Jan_28_2017
#macs_homer_overlap.sh

#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/macs_homer_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/macs_homer_e
#PBS -N macs_homer
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

cd /scratch.global/nosha003/schmitz/atacseq

# make file of all overlapping peaks (in MACs and HOMER peak calling)
perl /home/springer/nosha003/schmitz/scripts/macs_homer_overlap.pl -i /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}.pks.bed -I /scratch.global/nosha003/schmitz/atacseq/macs/${ID}_bowtie_macs_peaks.bed -o /scratch.global/nosha003/schmitz/atacseq/${ID}_macs_homer.txt

# count file
echo "homer"
wc -l /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}.pks.bed
echo "macs"
wc -l /scratch.global/nosha003/schmitz/atacseq/macs/${ID}_bowtie_macs_peaks.bed
echo "overlap"
wc -l /scratch.global/nosha003/schmitz/atacseq/${ID}_macs_homer.txt

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/macs_homer_overlap.sh