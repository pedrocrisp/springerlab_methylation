#Jan_28_2017
#macs_bins_intersect
#find overlap between macs/homer peaks and 100bp bins
	
#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=10gb
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

#sed '1d' counts/${ID}.bowtie.100bpcounts.bed > counts/${ID}.bowtie.100bpcounts.noheader.bed
sed '1d' /home/springer/nosha003/schmitz/atacseq/counts/${ID}_uniq.bowtie.100bpcounts.bed > counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed

#find intersect between macs peaks and 100bp bins
#bedtools intersect -a counts/${ID}.bowtie.100bpcounts.noheader.bed -b macs/${ID}_bowtie_macs_peaks.bed -bed -wa -wb > counts/${ID}_bowtie_macs_bins.bed
#bedtools intersect -a counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed -b macs/${ID}_bowtie_uniq_macs_peaks.bed -bed -wa -wb > counts/${ID}_uniq_bowtie_macs_bins.bed

#find intersect between homer peaks and 100bp bins
#bedtools intersect -a counts/${ID}.bowtie.100bpcounts.noheader.bed -b homer/peaks/${ID}.pks.bed -bed -wa -wb > counts/${ID}_bowtie_homer_bins.bed
bedtools intersect -b counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed -a /home/springer/nosha003/schmitz/atacseq/homer/${ID}_uniq_cutoff.pks.bed -bed -wa -wb > counts/${ID}_uniq_bowtie_homer_bins.bed

#bedtools intersect -b counts/B73_Leaf_rep1_uniq.bowtie.100bpcounts.noheader.bed -a /home/springer/nosha003/schmitz/atacseq/homer/B73_Leaf_rep1_uniq_cutoff.pks.bed -bed -wa -wb > counts/B73_Leaf_rep1_uniq_bowtie_homer_bins.bed
  
#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/peaks_bins_intersect.sh