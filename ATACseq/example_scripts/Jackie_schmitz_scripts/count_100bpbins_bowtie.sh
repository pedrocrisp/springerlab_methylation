#28Jan2017
#atacseq_count
#count reads for 100 bp tiles with bedtools intersect
	
#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_count_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_count_e
#PBS -N atacseq_count
#PBS -V
#PBS -r n

ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

#load modules
module load bedtools/2.25.0
module load samtools/1.3
	
cd /scratch.global/nosha003/schmitz/atacseq

#count number of reads aligning to each 100 bp window
#bedtools intersect -a align/${ID}_bowtie.clean.bam -b /home/springer/nosha003/database/B73v4_100bp_bins.gff -bed -wa -wb > counts/${ID}.bowtie.100bpbins.bed
bedtools intersect -a align/${ID}_bowtie_uniq.clean.bam -b /home/springer/nosha003/database/B73v4_100bp_bins.gff -bed -wa -wb > counts/${ID}_uniq.bowtie.100bpbins.bed
  
#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/count_100bpbins_bowtie.sh