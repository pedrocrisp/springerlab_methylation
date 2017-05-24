#!/bin/bash -l
#clean-up script
#Peter Crisp
#2017-05-24

# to run (USE CARE UNTESTED)
# cd <project folder>
# bash /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/09-clean-up.sh

# NOTE THIS FILES NEEDS TO BE UPDATED TO THE CURRENT PIPELINE

#removed trimmed fastq files
rm -rv analysis/trimmed/*.fq

#remove methratio.txt files because we have parsed these files to make BSMAP_out.txt
rm -rv analysis/BSMAPratio/*methratio.txt

# remove initial mapping files
# note these could be useful for looking at capture efficiency and distribution before filtering reads
rm -rv analysis/bsmapped

# remove intermediate mapping files
rm -rv analysis/bsmapped_filtered/*_sorted_MarkDup.bam
rm -rv analysis/bsmapped_filtered/*_sorted_MarkDup_pairs.bam
