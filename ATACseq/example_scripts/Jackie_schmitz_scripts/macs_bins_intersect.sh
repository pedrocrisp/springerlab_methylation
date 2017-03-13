#Jan_2_2017
#macs_bins_intersect
#find overlap between macs peaks and 100bp bins
	
#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/macs_bins_intersect_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/macs_bins_intersect_e
#PBS -N macs_bins_intersect
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

#load modules
module load bedtools/2.25.0

cd /scratch.global/nosha003/schmitz/atacseq

#find intersect between macs peaks and 100bp bins
bedtools intersect -a counts/${ID}.100bpcounts.bed -b macs/${ID}_macs_peaks.bed -bed -wa -wb > counts/${ID}_macs_bins_intersect.bed
  
#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/macs_bins_intersect.sh