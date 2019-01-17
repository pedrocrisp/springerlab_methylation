#!/bin/bash
#set -xe
set -xeuo pipefail

usage="USAGE:
bash 06-hotspot2-peaks_qsub.sh <sample_list.txt> <bam_folder> <chrom_sizes> <centers.starch>"

#define stepo in the pipeline - should be the same name as the script
step=06-hotspot2-peaks

######### Setup ################
sample_list=$1
bam_folder=$2
chrom_sizes=$3
centers_starch=$4

if [ "$#" -lt "4" ]
then
echo $usage
exit -1
else
echo "initiating peak calling jobs on $bam_folder folder"
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
-v LIST=${sample_list},bam_folder=${bam_folder},chrom_sizes=${chrom_sizes},centers_starch=${centers_starch} \
$script_to_qsub

# to run
# bash \
# /home/springer/pcrisp/gitrepos/springerlab_methylation/ATAC_seq/6-hotspot2-peaks_qsub.sh \
# samples.txt \
# align_bowtie2_k100_no_mismatch \
# ~/ws/refseqs/barley/hotspot/hotspot_barley_morex_IBSC_PGSB_v2_split_inc_Pt.chrom.sizes \
# ~/ws/refseqs/barley/hotspot/barley_morex_IBSC_PGSB_v2_split_inc_Pt.centers.starch
