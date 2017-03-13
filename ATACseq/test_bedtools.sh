#!/bin/bash
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /scratch.global/pcrisp/ATACseq_Jackie/methylation_ATACpeaks_bedtest_o
#PBS -e /scratch.global/pcrisp/ATACseq_Jackie/methylation_ATACpeaks_bedtest_e
#PBS -N methylation_ATACpeaks_bedtest
#PBS -V
#PBS -r n

#load modules
module load bedtools/2.25.0

#grab the relevant columns with cut and use sed to exclude the header
cut -f 1-3,4,8,12 /home/springer/pcrisp/ws/analysis/ATACseq_Jackie/data_Jackie/B73.tab  | sed '1d' > /scratch.global/pcrisp/ATACseq_Jackie/B73_tab.bed

# closest homer peak to each methylation bin
bedtools closest -D a -t first -a /scratch.global/pcrisp/ATACseq_Jackie/B73_tab.bed -b /home/springer/pcrisp/ws/analysis/ATACseq_Jackie/data_Jackie/B73_Leaf_uniq.pks.bed > /scratch.global/pcrisp/ATACseq_Jackie/homer_methyl_tabfile.txt

#make file for metaplot
perl /home/springer/nosha003/schmitz/scripts/methyl_homer_metaplot.pl -i /scratch.global/pcrisp/ATACseq_Jackie/homer_methyl_tabfile.txt -o /scratch.global/pcrisp/ATACseq_Jackie/homer_methyl_tabfile_metaplot.txt

#qsub /home/springer/pcrisp/ws/analysis/ATACseq_Jackie/scripts/test_bedtools.sh

