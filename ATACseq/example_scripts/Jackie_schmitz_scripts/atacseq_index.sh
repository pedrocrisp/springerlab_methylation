#25Jan2017
#atacseq_index.sh
#!/bin/bash

#Bowtie_Alignment
#PBS -l walltime=10:00:00,nodes=1:ppn=6,mem=10gb	
#PBS -o /scratch.global/nosha003/schmitz/e_o/atacseq_index_o
#PBS -e /scratch.global/nosha003/schmitz/e_o/atacseq_index_e
#PBS -N atacseq_index
#PBS -V
#PBS -r n

module load bowtie/0.12.7

bowtie-build /home/springer/nosha003/database/B73v4/Zea_mays.AGPv4.dna.genome.fa  /home/springer/nosha003/database/bowtie1db/Zea_mays.AGPv4.dna.genome

# qsub /home/springer/nosha003/schmitz/scripts/atacseq_index.sh