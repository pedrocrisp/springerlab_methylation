#28Jan2017
#macs_peaks
#count reads for 100 bp tiles with macs peak finding algorithm
	
#!/bin/bash
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/nosha003/schmitz/atacseq/e_o/macs_o
#PBS -e /scratch.global/nosha003/schmitz/atacseq/e_o/macs_e
#PBS -N macs
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"


#load modules
module load macs/1.4.1
module load R

cd /scratch.global/nosha003/schmitz/atacseq

#find peaks with macs peak finding algorithm
#macs14 callpeak -t align/${ID}_bowtie.clean.bam -f BAM -n macs/${ID}_bowtie_macs

# uniq
macs14 callpeak -t align/${ID}_bowtie_uniq.clean.bam -f BAM -n macs/${ID}_bowtie_uniq_macs
  
#qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/macs_bowtie.sh