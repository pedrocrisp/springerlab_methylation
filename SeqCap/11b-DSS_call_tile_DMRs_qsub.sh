#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 11b-DSS_call_tile_DMRs_qsub.sh <contrasts_list.txt> <DMR_contrasts_table_file> <path_to_data_files>
for example:
bash \
/home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/11b-DSS_call_tile_DMRs_qsub.sh \
DMR_tests_combos_all_list.txt \
DMR_tests_combos_all_table.tsv \
analysis/tiles_filtered_4C_2x
"

#define stepo in the pipeline - should be the same name as the script
step=11a-DSS_call_tile_DMRs

######### Setup ################
sample_list=$1
DMR_contrasts_table_file=$2
path_to_data_files=$3

if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "Submitting samples listed in '$sample_list' for filtering"
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
-v LIST=${sample_list},DMR_contrasts_table_file=$DMR_contrasts_table_file,path_to_data_files=$path_to_data_files \
$script_to_qsub
