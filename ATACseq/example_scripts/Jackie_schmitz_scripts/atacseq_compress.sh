#26_Jan_2017
#atacseq_compress.sh
#!/bin/bash

#PBS -l walltime=02:00:00,nodes=1:ppn=1,mem=2gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_compress_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_compress_e
#PBS -N atacseq_compress
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

#mkdir /home/springer/nosha003/schmitz/atacseq/fastq/compress

gzip /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE}.fastq > /home/springer/nosha003/schmitz/atacseq/fastq/compress/${ID}_1.fastq.gz

gzip /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE2}.fastq > /home/springer/nosha003/schmitz/atacseq/fastq/compress/${ID}_2.fastq.gz

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_compress.sh