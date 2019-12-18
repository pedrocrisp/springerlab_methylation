#!/bin/bash -l
#09-clean-up script
#Peter Crisp
#2017-05-24
# script to purge non-essential files to save space post running the pipeline
# copy files that should be backed up to desired dir

#############
# It could take some time so best to run in an interactive session in a screen eg
# screen -S copy_data
# qsub -I -l walltime=4:00:00,nodes=1:ppn=1,mem=4gb
# then to run: (USE CARE UNTESTED)
# cd <project folder> but can run from anywhere
# bash \
# /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/09-clean-up.sh \
# <full path to projectFolder>
# <"purge" or "no_purge"> \
# <"copy_to_home" or "no_copy_to_home">
# <full path to destination dir on home> \
# <"s3" or "no_s3"> \
# <destination on s3>

# # eg purge first:
# bash \
# /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/09-clean-up.sh \
# /scratch.global/pcrisp/m162_seqcap \
# purge \
# no_copy_to_home \
# NA \
# no_s3 \
# NA

# # inspect then copy to home and s3
# bash \
# /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/09-clean-up.sh \
# /scratch.global/pcrisp/m162_seqcap \
# no_purge \
# copy_to_home \
# ~/ws/analysis \
# s3 \
# s3://springer-pcrisp-WGBS/

#############

# The script can be run after step(s) 05 summarise methylation.
# NOTE step 05.2 context means needds to be run before this step

projectFolder=$1
purge=$2
copy_to_home=$3
home_dir_destination=$4
copy_to_s3=$5
s3_bucket=$6

cd $projectFolder

echo "Current folder is $PWD it will be $purge"
echo "Home copy? $copy_to_home"
echo "Home destionation $home_dir_destination"
echo "s3 copy? $copy_to_s3"
echo "s3 destination $s3_bucket"
echo "folder for copy to s3 $projectFolder"


# purge
if [ "$purge" == "purge" ]
then

# lets hope you backed up the reads somewhere!
rm -rv reads

#removed trimmed fastq files
rm -rv analysis/trimmed/*.fq

# remove initial mapping files
# note these could be useful for looking at capture efficiency and distribution before filtering reads
rm -rv analysis/bsmapped

# remove intermediate mapping files
rm -rv analysis/bsmapped_filtered/*_sorted_MarkDup.bam
rm -rv analysis/bsmapped_filtered/*_sorted_MarkDup_pairs.bam

#move HsMetrics and MarkDuplicates output log files
mkdir -p analysis/HsMetrics_deDups_logs
mv analysis/bsmapped_filtered/*.txt analysis/HsMetrics_deDups_logs/

###### BSMAPratio
# move bigWigs out of the BSMAPratio folder so they can be copied separately
mkdir -p analysis/BSMAPratio_bigWigs
mv -v analysis/BSMAPratio/*.bigWig analysis/BSMAPratio_bigWigs/
# it is recommended that the bigwigs are coppied to s3 and deleted from home... too big otherwise, they can be regenerated if required

# remaining files in BSMAPratio are quite large,
# it seems to make more sense to delete this folder now and keep the bams (half the size of *BSMAP_out.txt)
# ie these files could be easily recreated from bams
# *BSMAP_out.bg
# *BSMAP_out_subcontext.txt
# *BSMAP_out.txt
# *methratio.txt
rm -rv analysis/BSMAPratio

###### tiles
# I make tile bigwigs and also tile bed files, but I dont really use these
# the bed files may be usefull but they need coverage filtering first anyway
# so delete both BUT IF I DO WANT THEM OMIT THIS OR RERUN THE 07-tiles step
# the bg file is intermediate for making the bigWigs, it might also be useful for coverage analysis in the future FYI
rm -rv analysis/tiles/*.bg
rm -rv analysis/tiles/*.bigWig
rm -rv analysis/tiles/*.bed

else

  echo "Skiping purge"
fi

  if [ "$copy_to_home" == "copy_to_home" ]
  then

##################
# copy the remainder to home?
rsync -rhivPt $projectFolder $home_dir_destination/

# #sync log files
# rsync -rhivPt logs $home_dir_destination/
# #extra logs...
# rsync -rhivPt analysis/logs $home_dir_destination/
# #
# rsync -rhivPt analysis/ConversionRate $home_dir_destination/
# rsync -rhivPt analysis/fastqc $home_dir_destination/
# rsync -rhivPt analysis/mC_bigWigs $home_dir_destination/
# # BSMAPratio_bigWigs folder has to be created in the purge step above
# rsync -rhivPt analysis/BSMAPratio_bigWigs $home_dir_destination/
# rsync -rhivPt analysis/OnTargetCoverage $home_dir_destination/
# rsync -rhivPt analysis/tiles $home_dir_destination/
# rsync -rhivPt analysis/tiles_bigWigs $home_dir_destination/
# rsync -rhivPt analysis/analysis/HsMetrics_deDups_logs $home_dir_destination/
#
# #catch any auxillary text files such as sample list
# rsync -rhivPt *.txt $home_dir_destination/
#
# # BSMAPratio? this is the biggest folder! keep on s3...
# # rsync -rhivPt analysis/BSMAPratio $home_dir_destination/

else

  echo "Skiping copy to home"
fi

  if [ "$copy_to_s3" == "copy_to_s3" ]
  then

    ##################

    # THINGS THAT SHOULD BE KEPT AT LEAST ON S3 that will not be copied to home below:

    # The "BSMAPratio" folder will be large but that is a key output, save on s3 or on home
    # The "reads" folder should be redundant, but why not push a copy to s3 for backup...
    ## reads could be deleted if space is an issue
    # THE "bsmapped_filtered" folder has filtered bams that are used for all downstream analysis,
    ## "bsmapped_filtered" useful for IGV, it also has
    ## "bsmapped_filtered" also could be used to re-extract the methylation counts if required

    # it is recommended that everything left after the purge above should be copied to s3
    # eg
    # s3cmd sync --verbose SeqCap_2_McGinnis/ s3://springer-pcrisp-seqcap/SeqCap_2_McGinnis/
    # 76 GB took just under 2 hr to sync ultimately


    #get md5sum for each files
    find $projectFolder -type f -exec md5sum {} \; > md5sum.txt

    # sync
    s3cmd sync --verbose $projectFolder $s3_bucket

else

  echo "Skiping copy to s3"

  fi
