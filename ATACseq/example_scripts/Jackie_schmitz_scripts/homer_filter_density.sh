#9Feb2017
#homer_filter_density.sh
#!/bin/bash

#filter peaks
#PBS -l walltime=08:00:00,nodes=1:ppn=1,mem=2gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_filter_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_filter_e
#PBS -N atacseq_filter
#PBS -V
#PBS -r n

ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
READS="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 2)"

perl /home/springer/nosha003/schmitz/scripts/homer_filter_density.pl -i /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_uniq.pks.bed -o /scratch.global/nosha003/schmitz/atacseq/homer/peaks/${ID}_uniq_filter.pks.bed -p ${READS} -h help
#perl /home/springer/nosha003/schmitz/scripts/homer_filter_density.pl -i /scratch.global/nosha003/schmitz/atacseq/homer/peaks/B73_Leaf_rep1_uniq.pks.bed -o /scratch.global/nosha003/schmitz/atacseq/homer/peaks/B73_Leaf_rep1_uniq_filter.pks.bed -p 12197384 -h help

# qsub -t 1-17 -v LIST=/home/springer/nosha003/schmitz/scripts/atacseq_names_reads.txt /home/springer/nosha003/schmitz/scripts/homer_filter_density.sh