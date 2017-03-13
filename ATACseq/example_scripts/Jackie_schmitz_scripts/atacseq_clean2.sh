#24_Jan_2017
#atacseq_clean2.sh
#!/bin/bash

#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=2gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_clean_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_clean_e
#PBS -N atacseq_clean
#PBS -V
#PBS -r n

#sequencing data: /home/springer/nosha003/schmitz/atacseq/fastq

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load cutadapt/1.8.1
module load fastqc/0.10.1
module load python
module load fastx_toolkit
module load trimmomatic

#cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A ${ADAPTER} -o /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_1.fastq -p /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_2.fastq /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE}.fastq /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE2}.fastq
#fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_1.fastq
#fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_2.fastq

#cutadapt -a ${ADAPTER} -A ${ADAPTER} -o /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe_1.fastq -p /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe_2.fastq /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE}.fastq /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE2}.fastq
#fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe_1.fastq
#fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe_2.fastq

#try running cutadapt in SE mode instead
#cutadapt -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b ${ADAPTER} -f fastq -m 30 -q 10 --quality-base=33 -o /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe1.fastq /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE}.fastq
#cutadapt -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b ${ADAPTER} -f fastq -m 30 -q 10 --quality-base=33 -o /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe2.fastq /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE2}.fastq

#fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe1.fastq
#fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/clean/${ID}_pe2.fastq

#try trimmomatix 
java -Xmx7000M -jar $TRIMMOMATIC/trimmomatic.jar PE -basein /home/springer/nosha003/schmitz/atacseq/fastq/${SAMPLE}.fastq.gz -baseout /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}.fastq ILLUMINACLIP:$TRIMMOMATIC/adapters/all_illumina_adapters.fa:2:30:10:2:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35

fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_1P.fastq
fastqc -o /scratch.global/nosha003/schmitz/atacseq/fastqc --noextract -f fastq /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_2P.fastq

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_clean2.sh