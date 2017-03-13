#Jan_9_2017
#tdf_rna.sh
#!/bin/bash

#Bedgraph 
#PBS -l walltime=08:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/tdf_rna_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/tdf_rna_e
#PBS -N tdf_rna
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

module load igv
module load java

cd /scratch.global/nosha003/schmitz/rnaseq/align/tophat2

/home/springer/nosha003/igv/IGVTools/igvtools toTDF ${ID}.igv ${ID}.tdf /home/springer/nosha003/database/B73v4/Zea_mays.AGPv4.dna.genome.fa

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/rnaseq_names.txt /home/springer/nosha003/schmitz/scripts/tdf_rna.sh