#Jan_5_2017
#rnaseq_bedgraph.sh
#!/bin/bash

#Bedgraph 
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/bedgraph_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/bedgraph_e
#PBS -N bedgraph_rnaseq
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load bedtools

cd /scratch.global/nosha003/schmitz/rnaseq/align/tophat2
    
bedtools genomecov -ibam ${ID}.tophat2-g20.bam -bg > ${ID}.bedgraph

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/rnaseq_names.txt /home/springer/nosha003/schmitz/scripts/rnaseq_bedgraph.sh