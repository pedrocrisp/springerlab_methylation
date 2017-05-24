#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 03-fix-sort-bsmap_qsub.sh <sample_list.txt>
for example:
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/03-fix-sort-bsmap_qsub.sh \
single_sample.txt
"

#define stepo in the pipeline - should be the same name as the script
step=03-fix-sort-bsmap

######### Setup ################
sample_list=$1

if [ "$#" -lt "1" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for fixing and sorting"
cat $sample_list
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

########## Run #################

#make log and analysis folders
#make logs folder if it doesnt exist yet
mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)

#make analysis dir if it doesnt exist yet
analysis_dir=analysis
mkdir -p $analysis_dir

#make trimmgalore logs folder, timestamped
log_folder=logs/${timestamp}_${step}
mkdir $log_folder

#script path and cat a record of what was run
script_to_qsub=${scriptdir}/${step}.sh
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
qsub -t $qsub_t \
-o ${log_folder}/${step}_o \
-e ${log_folder}/${step}_e \
-v LIST=${sample_list} \
$script_to_qsub
