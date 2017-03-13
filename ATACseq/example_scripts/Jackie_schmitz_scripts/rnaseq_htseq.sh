#Jan_3_2017
#htseq.sh
#!/bin/bash

#Samstat_metrics
#PBS -l walltime=01:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/htseq_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/htseq_e
#PBS -N htseq
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

cd /scratch.global/nosha003/schmitz/rnaseq

#load modules
module load htseq/0.5.3
module load samtools/1.3
module load python/2.7.1
module load bedtools

#sort bam files
samtools sort -O bam -T temp align/tophat2/${ID}.tophat2-g20.bam > align/tophat2/${ID}.sort.bam
samtools index align/${ID}.sort.bam

#count number of reads aligning to each gene
#samtools view... change bam to sam file
#python to run HTseq
samtools view align/tophat2/${ID}.tophat2-g20.bam | python -m HTSeq.scripts.count -s no -t gene -i ID -m union -a 20 - /home/springer/nosha003/database/Zea_mays.AGPv4.32.gene.sort.gff3 > align/htseq/${ID}.htseqcount

#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/rnaseq_names.txt /home/springer/nosha003/schmitz/scripts/rnaseq_htseq.sh
#exel file of counts (align_summary)