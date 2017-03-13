#November_30_2016
#rnaseq_align.sh
#!/bin/bash

#TopHat2_Alignment
#PBS -l walltime=10:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/rnaseq_align_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/rnaseq_align_e
#PBS -N rnaseq_align
#PBS -V
#PBS -r n

#sequencing data: /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_047

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

if [ ! -d ${ID} ];
then
    mkdir ${ID}
fi

module load bowtie2/2.2.4
module load samtools/1.3
module load tophat/2.0.13

cd /scratch.global/nosha003/schmitz/rnaseq/clean
    
#tophat2 -p 6 -g 20 -i 5 -I 20000 -o /scratch.global/nosha003/schmitz/rnaseq/align/${ID}_pe_tophat2 /home/springer/nosha003/B73_RefGen_v4/v4_analysis/index/B73v4_pseudomolecules /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_047/${SAMPLE}.fastq /home/springer/data_release/umgc/hiseq/161110_D00635_0177_ACA37PANXX/Springer_Project_047/${SAMPLE2}.fastq

cp /scratch.global/nosha003/schmitz/rnaseq/align/${ID}_pe_tophat2/accepted_hits.bam /scratch.global/nosha003/schmitz/rnaseq/align/tophat2/${ID}.tophat2-g20.bam

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/rnaseq_names.txt /home/springer/nosha003/schmitz/scripts/rnaseq_align.sh