#Jan_9_2017
#count_igv.sh
#!/bin/bash

#Bedgraph 
#PBS -l walltime=08:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/igv_atac_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/igv_atac_e
#PBS -N igv_atacseq
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

cd /scratch.global/nosha003/schmitz/atacseq/align
    
perl /home/springer/nosha003/schmitz/scripts/bedgraph_igv.pl -i ${ID}.bedgraph -o ${ID}.igv

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/count_igv.sh