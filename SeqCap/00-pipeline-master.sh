#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 04-summarise_methylation_qsub.sh <sample_list.txt> <genome.fa> <CalculateHsMetrics_reference.bed> <intersect_regions.bed>
for example:
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/00-pipeline-master.sh \
single_sample.txt \
/home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa \
/home/springer/pcrisp/ws/refseqs/maize/seqcapv2_onTarget-for-picard.bed \
/home/springer/pcrisp/ws/refseqs/maize/BSseqcapv2_specific_regions.bed
"

#define stepo in the pipeline - should be the same name as the script
step=00-pipeline-master

######### Setup ################
sample_list=$1
genome_reference=$2
CalculateHsMetrics_reference=$3
intersect_regions_ref=$4
if [ "$#" -lt "4" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for SeqCap analysis"
cat $sample_list
echo genome reference is $genome_reference
echo CalculateHsMetrics_reference is $CalculateHsMetrics_reference
echo intersect regions are $intersect_regions_ref
fi

#number of samples
number_of_samples=`wc -l $sample_list | awk '{print $1}'`
if [[ "$number_of_samples" -eq 1 ]]
then
qsub_t=1
else
qsub_t="1-${number_of_samples}"
fi
echo "argument to be passed to qsub -t is '$qsub_t'"

#find script to run, makes it file system agnostic
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi

########## some folders and logs  #################

#make log and analysis folders
#make logs folder if it doesnt exist yet
mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)

#make analysis dir if it doesnt exist yet
analysis_dir=analysis
mkdir -p $analysis_dir

#script path and cat a record of what was run
log_folder0=logs/${timestamp}_00-pipeline-master
mkdir $log_folder0
cat $0 > ${log_folder0}/qsub_runner.log

########## scripts #################
step1=01-trim_galore
script_to_qsub1=${scriptdir}/${step1}.sh
log_folder1=logs/${timestamp}_${step1}
mkdir $log_folder1
cat $script_to_qsub1 > ${log_folder1}/script.log

step2=02-bsmap
script_to_qsub2=${scriptdir}/${step2}.sh
log_folder2=logs/${timestamp}_${step2}
mkdir $log_folder2
cat $script_to_qsub2 > ${log_folder2}/script.log

step3=03-filter
script_to_qsub3=${scriptdir}/${step3}.sh
log_folder3=logs/${timestamp}_${step3}
mkdir $log_folder3
cat $script_to_qsub3 > ${log_folder3}/script.log

step4=04-summarise_methylation
script_to_qsub4=${scriptdir}/${step4}.sh
log_folder4=logs/${timestamp}_${step4}
mkdir $log_folder4
cat $script_to_qsub4 > ${log_folder4}/script.log

########## qsub #################

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
FIRST=$(qsub -t $qsub_t \
-o ${log_folder1}/${step1}_o \
-e ${log_folder1}/${step1}_e \
-v LIST=${sample_list} \
$script_to_qsub1)
echo $FIRST

SECOND=$(qsub -W depend=afterany:$FIRST $script_to_qsub1 \
-t $qsub_t \
-o ${log_folder2}/${step2}_o \
-e ${log_folder2}/${step2}_e \
-v LIST=${sample_list},genome_reference=$genome_reference \
$script_to_qsub2)
echo $SECOND

THIRD=$(qsub -W depend=afterany:$SECOND $script_to_qsub2 \
-t $qsub_t \
-o ${log_folder3}/${step3}_o \
-e ${log_folder3}/${step3}_e \
-v LIST=${sample_list},genome_reference=$genome_reference,CalculateHsMetrics_reference=$CalculateHsMetrics_reference \
$script_to_qsub3)
echo $THIRD

FORTH=$(qsub -W depend=afterany:$THIRD $script_to_qsub3 \
-t $qsub_t \
-o ${log_folder4}/${step4}_o \
-e ${log_folder4}/${step4}_e \
-v LIST=${sample_list},genome_reference=$genome_reference,intersect_regions_ref=${intersect_regions_ref \
$script_to_qsub4)
echo $FORTH

##############################
