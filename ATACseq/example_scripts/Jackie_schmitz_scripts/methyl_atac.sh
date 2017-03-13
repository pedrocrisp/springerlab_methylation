#Jan_9_2017
#methyl_atac.sh
#!/bin/bash

#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/methyl_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/methyl_e
#PBS -N methyl
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
SAMPLE2="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"
ADAPTER="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 4)"

perl /home/springer/nosha003/schmitz/scripts/methylation_atacseq.pl -i /scratch.global/nosha003/schmitz/atacseq/counts/${ID}_macs_bins_counts2.txt -I /home/springer/nosha003/methylation_B73v4/B73_nonzero.tab -o /scratch.global/nosha003/schmitz/atacseq/methylation/${ID}_methyl_nonzero.txt

# qsub -t 1-12 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names.txt /home/springer/nosha003/schmitz/scripts/methyl_atac.sh
