#Dec_2_2016
#atacseq_count
#count reads for 100 bp tiles with bedtools intersect
	
#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/atacseq/e_o/atacseq_count_o
#PBS -e /scratch.global/nosha003/schmitz/atacseq/e_o/atacseq_count_e
#PBS -N atacseq_count
#PBS -V
#PBS -r n

ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

#load modules
module load bedtools
module load samtools
	
cd /scratch.global/nosha003/schmitz/atacseq

#sort and index alignment file
samtools sort -O bam -T tmp  align/${ID}_pe_bowtie2.bam > align/${ID}_sort.bam
samtools index align/${ID}_sort.bam

samtools sort -O sam -T tmp align/${ID}_pe_bowtie2.sam > align/${ID}_sort.sam

#count number of reads aligning to each 100 bp window
bedtools intersect -a align/${ID}_pe_bowtie2.bam -b /home/springer/nosha003/database/B73v4_100bp_bins.gff -bed -wa -wb > counts/${ID}.100bpbins.bed

## all bins
perl /home/springer/nosha003/wenli_data/scripts/bincounts.pl -i /home/springer/nosha003/database/B73v4_100bp_bins.gff -I counts/${ID}.100bpbins.bed -o counts/${ID}.counts.txt

## all bins bed count file
cut -f 13,16-21 counts/${ID}.100bpbins.bed | sort | uniq > counts/${ID}.100bpbins_cut.bed
perl /home/springer/nosha003/wenli_data/scripts/bincounts_bed.pl -i counts/${ID}.counts.txt -I counts/${ID}.100bpbins_cut.bed -o counts/${ID}.100bpcounts.bed

  
#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/100bp_counts.sh
