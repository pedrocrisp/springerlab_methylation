#30Nov2016

#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_fastqc_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_fastqc_e
#PBS -N atacseq_fastqc
#PBS -V
#PBS -r n

#sequencing data: /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_046

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

module load fastqc/0.10.1

fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_046/${SAMPLE}.fastq

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/atacseq_names.txt /home/springer/nosha003/schmitz/atacseq_fastqc.sh