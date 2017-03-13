#24Jan2017
#atacseq_align2.sh
#!/bin/bash

#Bowtie_Alignment
#PBS -l walltime=10:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_align2_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_align2_e
#PBS -N atacseq_align2
#PBS -V
#PBS -r n

#sequencing data: /home/springer/nosha003/schmitz/atacseq/fastq

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load bowtie/2.3.0
module load bowtie/1.1.2
module load samtools/1.3

cd /scratch.global/nosha003/schmitz/atacseq

# trimmomatic

## bowtie1
#bowtie /home/springer/nosha003/database/bowtie1db/Zea_mays.AGPv4.dna.genome -X 1000 -1 /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_1P.fastq -2 /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_2P.fastq -S /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_trimmomatic.sam -t -p 4 -v 0  --best --strata -m 3

bowtie /home/springer/nosha003/database/bowtie1db/Zea_mays.AGPv4.dna.genome -X 1000 -1 /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_1P.fastq -2 /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_2P.fastq -S /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie_trimmomatic_uniq.sam -t -p 4 -v 0 --best --strata -m 1

## bowtie2
#bowtie2 -p 6 --phred33 -x /home/springer/nosha003/B73_RefGen_v4/v4_analysis/index/B73v4_pseudomolecules -1 /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_1P.fastq -2 /scratch.global/nosha003/schmitz/atacseq/trimmomatic/${ID}_2P.fastq -S /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie2_trimmomatic.bam
samtools view -S -b -o /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie2_trimmomatic.sam /scratch.global/nosha003/schmitz/atacseq/align/${ID}_bowtie2_trimmomatic.bam
    # filter based on quality score to do same as bowtie1 -m 3 

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_align2.sh