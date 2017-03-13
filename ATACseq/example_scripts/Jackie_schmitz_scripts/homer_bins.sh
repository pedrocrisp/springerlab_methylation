#Feb_24_2017
#homer_bins
#find overlap between homer peaks and 100bp bins (only keep complete overlap of bins and peaks)
	
#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=10gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/homer_bins_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/homer_bins_e
#PBS -N homer_bins
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

#load modules
module load bedtools/2.25.0

# intersect 100bp counts bins with homer peaks and keep only bins completely overlapping with peaks
sed '1d' /home/springer/nosha003/schmitz/atacseq/counts/${ID}_uniq.bowtie.100bpcounts.bed > /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed

bedtools intersect -b /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed -a /home/springer/nosha003/schmitz/atacseq/homer/${ID}_uniq_cutoff.pks.bed -bed -wa -wb > /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq_bowtie_homer_bins.bed

perl /home/springer/nosha003/schmitz/scripts/bincounts_bed_homer.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq_bowtie_homer_bins.bed -o /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq_bowtie_homer_bins_counts_completeoverlap.txt

#scp nosha003@login.msi.umn.edu:/scratch.global/nosha003/schmitz/atacseq/counts/*_uniq_bowtie_homer_bins_counts_completeoverlap.txt /Volumes/CBS/Groups/LAB-springer/Jaclyn/Graduate_projects/schmitz_data/atacseq/methylation/.

# go through all 100bp bins and declare homer peak (YES/NO) in additional column
perl /home/springer/nosha003/schmitz/scripts/bincounts_homer_calls.pl -i /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq.bowtie.100bpcounts.noheader.bed -I /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_uniq_bowtie_homer_bins_counts_completeoverlap.txt -o /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_bincounts_homercall.txt

#scp nosha003@login.msi.umn.edu:/scratch.global/nosha003/schmitz/atacseq/counts/B73_Leaf_rep1_bincounts_homercall.txt /Volumes/CBS/Groups/LAB-springer/Jaclyn/Graduate_projects/schmitz_data/atacseq/methylation/.

#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/homer_bins.sh