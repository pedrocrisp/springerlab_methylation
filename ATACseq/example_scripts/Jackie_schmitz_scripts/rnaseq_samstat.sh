#Jan_3_2017
#samstat.sh
#!/bin/bash

#Samstat_metrics
#PBS -l walltime=08:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/samstat_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/samstat_e
#PBS -N samstat
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

cd /scratch.global/nosha003/schmitz/rnaseq/samstat

#samstat
module load samtools/0.1.18

home/springer/nosha003/sourceforge.net/projects/samstat/files/samstat-1.5.1/src/samstat /scratch.global/nosha003/schmitz/rnaseq/align/tophat2/${ID}.tophat2-g20.bam

#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/rnaseq_names.txt /home/springer/nosha003/schmitz/scripts/rnaseq_samstat.sh