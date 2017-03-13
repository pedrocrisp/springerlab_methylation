#Jan_28_2017
#matrix_fpkm

#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=20gb
#PBS -o /scratch.global/nosha003/schmitz/e_o/matrix_fpkm_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/matrix_fpkm_e
#PBS -N matrix_fpkm
#PBS -V
#PBS -r n

SAMPLE="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 1)"
ID="$(/bin/sed -n ${PBS_ARRAYID}p ${LIST} | cut -f 3)"

cd /scratch.global/nosha003/schmitz/atacseq

## uniq 
# create matrix of atacseq counts
perl /home/springer/nosha003/schmitz/scripts/atacseq_matrix_uniq.pl --bins /home/springer/nosha003/database/B73v4_100bp_bins.gff --matrix_out /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_uniq.txt
# combine reps and correct for read depth
perl /home/springer/nosha003/schmitz/scripts/atac_rep_combine_uniq.pl -i /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_uniq.txt -o /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq.txt

# uniq combined reps
perl /home/springer/nosha003/schmitz/scripts/atacseq_matrix_uniq_combinereps.pl --bins /home/springer/nosha003/database/B73v4_100bp_bins.gff --matrix_out /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_uniq_combinerep.txt
# combine reps and correct for read depth
perl /home/springer/nosha003/schmitz/scripts/atac_combinerep_fpkm_uniq.pl -i /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_uniq.txt -o /home/springer/nosha003/schmitz/atacseq/counts/atacseq_counts_matrix_fpkm_uniq_combinerep.txt

# qsub /home/springer/nosha003/schmitz/scripts/atac_fpkm_uni.sh