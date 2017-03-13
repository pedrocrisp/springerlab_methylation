#Jan_5_2017
#atacseq_bedgraph.sh
#!/bin/bash

#Bedgraph 
#PBS -l walltime=08:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/bedgraph_atac_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/bedgraph_atac_e
#PBS -N bedgraph_atacseq
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load bedtools

cd /scratch.global/nosha003/schmitz/atacseq/align
    
bedtools genomecov -ibam ${ID}_sort.bam -bg > ${ID}.bedgraph

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/atacseq_bedgraph.sh