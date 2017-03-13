#November_30_2016
#atacseq_align.sh
#!/bin/bash

#Bowtie2_Alignment
#PBS -l walltime=10:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_align_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_align_e
#PBS -N atacseq_align
#PBS -V
#PBS -r n

#sequencing data: /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_046

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load bowtie2/2.2.4
module load samtools/1.3

cd /scratch.global/nosha003/schmitz/atacseq

bowtie2 -p 6 --phred33 -x /home/springer/nosha003/B73_RefGen_v4/v4_analysis/index/B73v4_pseudomolecules -1 /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_046/${SAMPLE}.fastq -2 /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_046/${SAMPLE2}.fastq -S /scratch.global/nosha003/schmitz/atacseq/align/${ID}_pe_bowtie2

samtools view -S -b -o /scratch.global/nosha003/schmitz/atacseq/align/${ID}_pe_bowtie2.bam  /scratch.global/nosha003/schmitz/atacseq/align/${ID}_pe_bowtie2

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_align.sh

