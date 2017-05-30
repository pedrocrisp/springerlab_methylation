#!/bin/bash -l
#clean-up script
#Peter Crisp
#2017-05-24
# script to purge non-essential files to save space post running the pipeline
# copy files that should be backed up to desired dir

#############
# It could take some time so best to run in an interactive session in a screen eg
# screen -S copy_data
# qsub -I -l walltime=4:00:00,nodes=1:ppn=1,mem=4gb
# then to run: (USE CARE UNTESTED)
# cd <project folder>
# bash \
# /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/09-clean-up.sh
# <destination on home>
#############

# NOTE THIS FILE NEEDS TO BE UPDATED TO THE CURRENT PIPELINE

home_dir_destination=$1

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

#move HsMetrics and MarkDuplicates output log files
mkdir -p analysis/HsMetrics_deDups_logs
mv analysis/bsmapped_filtered/*.txt analysis/HsMetrics_deDups_logs/

##################

# THINGS THAT SHOULD BE KEPT AT LEAST ON S3 that will not be copied to home below:

# The "BSMAPratio" folder will be large but that is a key output, save on s3 or on home
# The "reads" folder should be redundant, but why not push a copy to s3 for backup...
## reads could be deleted if space is an issue
# THE "bsmapped_filtered" folder has bams - useful for IGV, it also has

# it is recommended that everything left after the purge above should be copied to s3
# eg
# s3cmd sync --verbose SeqCap_2_McGinnis/ s3://springer-pcrisp-seqcap/SeqCap_2_McGinnis/
# 76 GB took just under 2 hr to sync ultimately

##################
# copy the important stuff to backed up file system
# NOTE DECIDE WHAT IS IMPORTANT AND HOW TO SEPERATE?

#sync log files
rsync -rhivPt logs $home_dir_destination/
#extra logs...
rsync -rhivPt analysis/logs $home_dir_destination/
#
rsync -rhivPt analysis/ConversionRate $home_dir_destination/
rsync -rhivPt analysis/fastqc $home_dir_destination/
rsync -rhivPt analysis/mC_bigWigs $home_dir_destination/
rsync -rhivPt analysis/OnTargetCoverage $home_dir_destination/
rsync -rhivPt analysis/tiles $home_dir_destination/
rsync -rhivPt analysis/tiles_bigWigs $home_dir_destination/
rsync -rhivPt analysis/analysis/HsMetrics_deDups_logs $home_dir_destination/

#catch any auxillary text files such as sample list
rsync -rhivPt *.txt $home_dir_destination/

# BSMAPratio?
rsync -rhivPt analysis/BSMAPratio $home_dir_destination/
