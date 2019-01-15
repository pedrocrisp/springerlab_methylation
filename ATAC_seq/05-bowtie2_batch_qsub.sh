#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 05-bowtie2_batch_qsub.sh <sample_list.txt> <reads_folder> <bt2_threads> <bt2_genome.fa> <multimapping_rate>"
#eg /home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa

#define stepo in the pipeline - should be the same name as the script
step=05-bowtie2_batch

######### Setup ################
sample_list=$1
reads_folder=$2
bt2_threads=$3
bt2_genome=$4
multimapping_rate=$5
if [ "$#" -lt "5" ]
then
echo $usage
exit -1
else
echo "initiating bowtie jobs on $reads_folder folder, bowtie2 can use $bt2_threads threads"
cat $sample_list
echo genome reference is $bt2_genome
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
-v LIST=${sample_list},reads_folder=$reads_folder,bt2_threads=$bt2_threads,bt2_genome=$bt2_genome,multimapping_rate=$multimapping_rate \
$script_to_qsub

# to run
# bash /home/springer/pcrisp/gitrepos/NGS-pipelines/smallRNAseqPipe1/05-bowtie2_batch_qsub.sh <sample_list.txt> <reads_folder> <bt2_threads> <bt2_genome.fa> <multimapping_rate>
# eg
# bash /home/springer/pcrisp/gitrepos/springerlab_methylation/SeqCap/02-bsmap_qsub.sh samples.txt /home/springer/pcrisp/ws/refseqs/maize/Zea_mays.AGPv4.dna.toplevel.fa
