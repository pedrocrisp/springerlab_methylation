#24Jan2017
#atacseq_bam.sh
#!/bin/bash

#sam to bam, sort, remove clonal, bam to bed
#PBS -l walltime=02:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_bam_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_bam_e
#PBS -N atacseq_bam
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load bowtie2/2.2.4
module load samtools/1.3
module load bedtools

cd /scratch.global/nosha003/schmitz/atacseq

#samtools view -bS /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_trimmomatic.sam > /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.bam
                 
#samtools sort /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.bam > /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.sorted
                 
#samtools rmdup /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.sorted  /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.clean.bam

#bedtools bamtobed -i /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.clean.bam > /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie.bed


#uniq
samtools view -bS /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_trimmomatic_uniq.sam > /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.bam
                 
samtools sort /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.bam  > /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.sorted
                 
samtools rmdup /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.sorted /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.clean.bam

bedtools bamtobed -i /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.clean.bam > /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_uniq.bed

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_bam.sh