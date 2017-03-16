#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 01-trim_galore_qsub.sh <sample_list.txt>"

######### Setup ################
sample_list=$1
if [ "$#" -lt "1" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for trimming"
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

########## Run #################

#make log and analysis folders
#make logs folder if it doesnt exist yet
mkdir -p logs

timestamp=$(date +%Y%m%d-%H%M%S)

#make analysis dir if it doesnt exist yet
analysis_dir=analysis
mkdir -p $analysis_dir

#make trimmgalore logs folder, timestamped
log_folder=logs/logs_trim_galore_${timestamp}
mkdir $log_folder

#script path and cat a record of what was run
script_dir=/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap
script_to_qsub=${script_dir}/01-trim_galore.sh
cat $script_to_qsub > ${log_folder}/script.log
cat $0 > ${log_folder}/qsub_runner.log

#submit qsub and pass args
#-o and -e pass the file locations for std out/error
#-v additional variables to pass to the qsub script including the PBS_array list and the dir structures
qsub -t $qsub_t \
-o ${log_folder}/trim_galore_o \
-e ${log_folder}/trim_galore_e \
-v LIST=${sample_list} \
$script_to_qsub

#to run
#bash /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/01-trim_galore_qsub.sh <sample_list.txt>
