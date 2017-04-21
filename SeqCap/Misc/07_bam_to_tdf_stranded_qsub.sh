#!/bin/bash
#set -xe
set -xeuo pipefail

#draft qsub script for making bigWigs using old pipeline

usage="USAGE:
bash 07_bam_to_tdf_stranded_qsub.sh <alignment folder> <number of threads> <strandedness of library> <chromosome.sizes file> <path_to_07_bam_to_tdf_stranded.sh>"
#eg /home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa

#define step in the pipeline - should be the same name as the script
step=07_bam_to_tdf_stranded

######### Setup ################
alignment_folder=$1
threads=$2
strandedness=$3
chromosome_sizes=$4
07_bam_to_tdf_stranded=$5
if [ "$#" -lt "4" ]
then
echo $usage
exit -1
else
echo Converting bams in $1 to tdfs and stranded bigWigs
echo Parellel with will using $2 threads
echo strandedness = $3
echo chromosome_sizes file $4
fi

#find script to run, makes it file system agnostic
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi

########## Run #################

#make log and analysis folders
#make logs folder if it doesnt exist yet
mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)

#make analysis dir if it doesnt exist yet
analysis_dir=analysis
mkdir -p $analysis_dir

#make tdf logs folder, timestamped
log_folder=logs/${timestamp}_${step}
mkdir ../$log_folder

#script path and cat a record of what was run
script_to_qsub=${07_bam_to_tdf_stranded}
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
qsub \
-o ${log_folder}/${step}_o \
-e ${log_folder}/${step}_e \
-v alignment_folder=${alignment_folder},genome_reference=$genome_reference \
$script_to_qsub
